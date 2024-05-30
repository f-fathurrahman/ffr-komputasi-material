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
gdml = GDMLPredict(model)

R, E_true, F_true = load_dataset(idx_data=10)
E_pred, F_pred = gdml.predict(R[None,:])

print("E_pred = ", E_pred)
print("E_true = ", E_true)
