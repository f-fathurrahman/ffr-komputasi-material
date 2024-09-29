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
print("^^^^ Initializing model parameters")
variables = model.init(
    random.PRNGKey(np.random.randint(42)),
    batch["coord"][0],
    batch["box"][0],
    static_args
)
print("^^^^ Done initializing model parameters")
# variables is NN parameters
print("Model initialized. Precision: %s. Parameter count: %d." % 
            ({"default": "fp32", "low": "fp32-16", "high": "fp64"}[INPUT.precision], 
            sum(i.size for i in tree_util.tree_flatten(variables)[0])))



# initialize optimizer
lr_scheduler = optax.exponential_decay(
    init_value=INPUT.lr,
    transition_steps=INPUT.decay_steps,
    decay_rate=(INPUT.lr_limit/INPUT.lr)**(INPUT.decay_steps/(INPUT.total_steps-INPUT.decay_steps)),
    transition_begin=0,
    staircase=True
)
optimizer = optax.adam(learning_rate=lr_scheduler, b2=INPUT.beta2)

print("Start optimizer init variables")
opt_state = optimizer.init(variables)
print("End optimizer init variables")

loss, loss_and_grad = model.get_loss_fn()
print("# Optimizer initialized, lr starts from %.1e. Starting training..." % INPUT.lr)

state = {"loss_avg":0., "iteration":0} | ({} if INPUT.model_type == "atomic" else {"le_avg":0., "lf_avg":0.})


@partial(jit, static_argnums=(4,))
def train_step(batch, variables, opt_state, state, static_args):
    print("*** ENTER train_step")
    r = lr_scheduler(state["iteration"]) / INPUT.lr
    if INPUT.model_type == "energy":
        pref = {
            "e": INPUT.s_pref_e*r + INPUT.l_pref_e*(1-r),
            "f": INPUT.s_pref_f*r + INPUT.l_pref_f*(1-r)
        }
        (loss_total, (loss_e, loss_f)), grads = loss_and_grad(variables, batch, pref, static_args)
        for key, value in zip(["loss_avg", "le_avg", "lf_avg"], [loss_total, loss_e, loss_f]):
            state[key] = state[key] * (1-1/INPUT.l_smoothing) + value
    else:
        loss_total, grads = loss_and_grad(variables, batch, static_args)
        state["loss_avg"] = state["loss_avg"] * (1-1/INPUT.l_smoothing) + loss_total
    updates, opt_state = optimizer.update(grads, opt_state, variables)
    variables = optax.apply_updates(variables, updates)
    state["iteration"] += 1
    print("*** EXIT train_step")
    return variables, opt_state, state

@partial(jit, static_argnums=(2,))
def val_step(batch, variables, static_args):
    if INPUT.model_type == "energy":
        pref = {"e": 1, "f": 1}
        _, (loss_e, loss_f) = loss(variables, batch, pref, static_args)
        return loss_e, loss_f
    else:
        loss_total = loss(variables, batch, static_args)
        return loss_total

tic = time()
for iteration in range(INPUT.total_steps+1):
    batch, type_count, lattice_args = train_data.get_batch(INPUT.label_bs, "label")
    static_args = nn.FrozenDict({"type_count":tuple(type_count), "lattice":lattice_args})
    variables, opt_state, state = train_step(batch, variables, opt_state, state, static_args)
    if iteration % INPUT.print_every == 0:
        if INPUT.use_val_data:
            val_batch, type_count, lattice_args = val_data.get_batch(INPUT.val_label_bs, "label")
            static_args = nn.FrozenDict({"type_count":tuple(type_count), "lattice":lattice_args})
            loss_val = val_step(val_batch, variables, static_args)
        beta = INPUT.l_smoothing * (1 - (1/INPUT.l_smoothing)**(iteration+1))
        print("Iter %7d" % iteration
              + " L %7.5f" % (state["loss_avg"]/beta)**0.5
              + (" LE %7.5f" % (state["le_avg"]/beta)**0.5 if INPUT.model_type == "energy" else "")
              + (" LF %7.5f" % (state["lf_avg"]/beta)**0.5 if INPUT.model_type == "energy" else "")
              + (" LEval %7.5f" % loss_val[0]**0.5 if INPUT.model_type == "energy" and INPUT.use_val_data else "")
              + (" LFval %7.5f" % loss_val[1]**0.5 if INPUT.model_type == "energy" and INPUT.use_val_data else "")
              + (" Lval %7.5f" % loss_val**0.5 if INPUT.model_type == "atomic" and INPUT.use_val_data else "")
              + " Time %.2fs" % (time()-tic))
        tic = time()
if INPUT.compress:
    model, variables = utils.compress_model(model, variables, INPUT.compress_Ngrids, INPUT.compress_rmin)
utils.save_model(INPUT.save_name, model, variables)
T = int(time() - TIC)
print("# Training finished in %dh %dm %ds." % (T//3600,(T%3600)//60,T%60))

