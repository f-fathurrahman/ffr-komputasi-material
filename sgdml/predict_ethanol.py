import numpy as np
from my_sgdml.predict import GDMLPredict
from my_sgdml.utils import io

model = np.load("m_ethanol.npz")
gdml = GDMLPredict(model)

#r,_ = io.read_xyz("ethanol01.xyz") # 9 atoms
r,_ = io.read_xyz("ethanol.xyz") # 9 atoms
e,f = gdml.predict(r)

print(r.shape) # (1,27)
print(e.shape) # (1,)
print(f.shape) # (1,27)


print(r.reshape(-1,3))
print(f.reshape(-1,3))
print(e)
