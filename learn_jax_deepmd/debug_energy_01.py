import INPUT01 as INPUT

#from INPUT01 import *

# override some input parameters
INPUT.total_steps = 1
INPUT.print_every = 1


# From here on you don"t need to change anything unless you know what you are doing
import numpy as np
from jax import jit, random, tree_util
import jax, optax, datetime
import flax.linen as nn
from time import time
from functools import partial

from my_deepmd_jax import data, utils
from my_deepmd_jax.dpmodel import DPModel


# From here on you don"t need to change anything unless you know what you are doing
import jax
import datetime
from time import time

from my_deepmd_jax import data

# This is not yet needed for loading data
if INPUT.precision == "default":
    print("Using default precision float32")
    jax.config.update("jax_default_matmul_precision", "float32")
if INPUT.precision == "high":
    print("Using float64")
    jax.config.update("jax_enable_x64", True)


TIC = time()
print("Program start at", datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"), "on device:", jax.devices()[:1])

if INPUT.model_type == "energy":
    labels = ["coord","box"] + ["force","energy"]
else:
    labels = ["coord","box"] + [INPUT.atomic_label]

print("labels = ", labels)

# Load data
train_data = data.DPDataset(INPUT.train_paths, labels, {"atomic_sel": INPUT.atomic_sel})

if INPUT.use_val_data:
    val_data = data.DPDataset(INPUT.val_paths, labels, {"atomic_sel": INPUT.atomic_sel})
else:
    val_data = None

# Compute lattice vectors?
train_data.compute_lattice_candidate(INPUT.rcut)
if INPUT.use_val_data:
    val_data.compute_lattice_candidate(INPUT.rcut)


params = {
    "embed_widths": INPUT.embed_widths[:-1] if INPUT.use_mp else INPUT.embed_widths,
    "embedMP_widths": INPUT.embed_widths[-1:] + INPUT.embedMP_widths if INPUT.use_mp else None,
    "fit_widths": INPUT.fit_widths,
    "axis": INPUT.axis_neurons,
    "Ebias": train_data.fit_energy() if INPUT.model_type == "energy" else None,
    "rcut": INPUT.rcut,
    "use_2nd": INPUT.use_2nd_tensor,
    "use_mp": INPUT.use_mp,
    "atomic":True if INPUT.model_type == "atomic" else False,
    "nsel": INPUT.atomic_sel if INPUT.model_type == "atomic" else None,
    "out_norm": 1. if INPUT.model_type == "energy" else train_data.get_atomic_label_scale()
}

data_stats = train_data.get_stats(INPUT.rcut, INPUT.getstat_bs)
model = DPModel(params|data_stats) # | will concat two dicts
print("Model params:", model.params)


# initialize model variables
batch, type_count, lattice_args = train_data.get_batch(1) # get one data point
static_args  = nn.FrozenDict({"type_count": type_count, "lattice": lattice_args})
variables = model.init(
    random.PRNGKey(np.random.randint(42)),
    batch["coord"][0],
    batch["box"][0],
    static_args
)


coord_N3 = batch["coord"][0]
box_33 = batch["box"][0]
nbrs_nm = None

# prepare input parameters
coord_N3, type_count, mask, compress, K, nsel, nbrs_nm = model.get_input(coord_N3, static_args, nbrs_nm)

A = model.params["axis"]
L = static_args["lattice"]["lattice_max"] if nbrs_nm is None else None

# compute relative coordinates x_3NM, distance r_NM, s(r) and normalized s(r)
x_n3m, r_nm = utils.get_relative_coord(coord_N3, box_33, type_count, static_args.get("lattice",None), nbrs_nm)
sr_nm = [ [utils.sr(r, model.params["rcut"]) for r in R] for R in r_nm ]
sr_norm_nm = [[r/std for r in R] for R,std in zip(sr_nm, model.params["sr_std"])]
sr_centernorm_nm = [[(r-mean)/std for r in R] for R,mean,std in zip(sr_nm, model.params["sr_mean"], model.params["sr_std"])]

# environment matrix: sr_norm_nm (0th-order), R_n3m (1st-order), R2_n6m (2nd-order)
x_norm_n3m = [
    [x/(r+1e-16)[:,None] for x,r in zip(X,R)] \
        for X,R in zip(x_n3m,r_nm)
]

R_n3m = [
    [3**0.5 * sr[:,None] * x for sr,x in zip(SR,X)] \
        for SR,X in zip(sr_norm_nm,x_norm_n3m)
]

R_n4m = [
    [utils.concat([sr[:,None],r], axis=1) for sr,r in zip(SR,R)] \
        for SR,R in zip(sr_norm_nm,R_n3m)
]

R_nsel6m = [
    [3*sr[:,None]*utils.tensor_3to6(x,axis=1,bias=1/3) for sr,x in zip(sr_norm_nm[nsel[i]], x_norm_n3m[nsel[i]])] \
        for i in range(len(nsel))
]

R_nselXm = [
    [ utils.concat([sr[:,None],r3] + ([r6] if model.params["use_2nd"] else []), axis=1) \
     for sr,r3,r6 in zip(sr_norm_nm[nsel[i]],R_n3m[nsel[i]],R_nsel6m[i])] \
        for i in range(len(nsel))
]


#Nembed = 0
#i = 0
#sri = sr_centernorm_nm[nsel[i]]
#rxi = R_nselXm[i]
#sr = sri[0]
#rx = rxi[0]
#net_name = "embedding_net_" + str(Nembed)
#net_vars = variables["params"][net_name]
#utils.embedding_net(model.params["embed_widths"]).apply({"params": net_vars}, sr[:,:,None], compress, rx)


Nembed = 0
list_embed = []
for i in range(len(nsel)):
    for sr, rx in zip(sr_centernorm_nm[nsel[i]], R_nselXm[i]):
        net_name = "embedding_net_" + str(Nembed)
        net_vars = variables["params"][net_name]
        list_embed.append(
            utils.embedding_net(model.params['embed_widths']).apply(
                {"params": net_vars}, sr[:,:,None], compress, rx
            )
        )
        Nembed += 1

T_NselXW = utils.concat(list_embed, K=K) / model.params['Nnbrs']


#T_NselW = T_NselXW[:,0] + model.param("Tbias", utils.zeros_init, T_NselXW.shape[-1:])
T_NselW = T_NselXW[:,0] + variables["params"]["Tbias"]
T_Nsel3W = T_NselXW[:,1:4]
T_Nsel6W = T_NselXW[:,4:]
G_NselAW = T_NselW[:,None]*T_NselW[:,:A,None] + (T_Nsel3W[:,:,None]*T_Nsel3W[:,:,:A,None]).sum(1)

# XXX only for use_2nd
G2_axis_Nsel6A = utils.tensor_3to6(T_Nsel3W[:,:,A:2*A], axis=1) + T_Nsel6W[:,:,A:2*A]
G_NselAW += (G2_axis_Nsel6A[...,None] * T_Nsel6W[:,:,None]).sum(1)

#fit_n1 = [utils.fitting_net(model.params["fit_widths"])(G) for G in utils.split(G_NselAW.reshape(G_NselAW.shape[0],-1),type_count,0,K=K)]
G_split = utils.split(G_NselAW.reshape(G_NselAW.shape[0],-1),type_count,0,K=K)
fit_n1 = []
Nfit = 0
for G in G_split:
    net_name = "fitting_net_" + str(Nfit)
    print("net_name = ", net_name)
    net_vars = variables["params"][net_name]
    print(G.shape)
    res = utils.fitting_net(model.params["fit_widths"]).apply({"params": net_vars}, G)
    print("res = ", res)
    #fit_n1.append(
    #    utils.fitting_net(model.params["fit_widths"]).apply({"params": net_vars}, G)
    #)
    Nfit += 1

#print("len of fit_n1 = ", len(fit_n1))
#print("len of mask = ", len(mask))
#print("K = ", K)
#
#list_E = []
#for f,Eb in zip(fit_n1, model.params["Ebias"]):
#    list_E.append(f[:,0] + Eb)
#print("list_E = ", list_E)

#pred = (mask * utils.concat(list_E, K=K)).sum()
#energy_pred = pred * model.params["out_norm"]


