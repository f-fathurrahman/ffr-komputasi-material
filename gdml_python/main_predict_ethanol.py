import numpy as np
#from my_sgdml.predict import GDMLPredict
from my_predict import GDMLPredict

def load_dataset(idx_data=0):
    dataset = np.load("DATASET/ethanol_dft.npz")
    R = dataset["R"][idx_data]
    E = dataset["E"][idx_data]
    F = dataset["F"][idx_data]
    return R, E, F

model = np.load("m_ethanol.npz")
#model = np.load("ethanol-aims.PBE.TS.light.tier.1-train200-sym6.npz")
gdml = GDMLPredict(model)

#idxs_train = model["idxs_train"]
R, E_true, F_true = load_dataset(idx_data=1)
#E_pred, F_pred = gdml.predict(R[None,:])
Natoms = R.shape[0]
print("Natoms = ", Natoms)
r = R.reshape(1,Natoms*3) # we need to reshape it
E_pred, F_pred = gdml.predict(r)

print("E_pred = ", E_pred)
print("E_true = ", E_true)
dE = abs(E_pred - E_true)
print("dE = ", dE)
print("Relative diff (in percent) = ", (100*dE/abs(E_true)))

Natoms = R.shape[0]
print("Natoms = ", Natoms)
F_pred = F_pred.reshape((Natoms,3))
print("F_pred = ")
print(F_pred)
print("F_true = ")
print(F_true)
dF = np.abs(F_pred - F_true)
print("dF = ")
print(dF)
