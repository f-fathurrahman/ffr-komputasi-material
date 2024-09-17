import torch
import numpy as np

coord = np.array(np.load("../DATASET/DeepPot/input_coord.npy", "r"), dtype="float32")
atom_types = np.loadtxt("../DATASET/DeepPot/type.raw", dtype=int)
elems = np.unique(atom_types).tolist()
atom_types = np.array([[elems.index(i) for i in atom_types]])
atom_types = atom_types.repeat(len(coord), axis=0)
grad = np.array(np.load("../DATASET/DeepPot/input_grad.npy", "r"), dtype="float32")

# Now, convert them to Tensor
#device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
device = torch.device("cpu")
coord = torch.from_numpy(coord).to(device)
atom_types = torch.from_numpy(atom_types).to(device)
grad = torch.from_numpy(grad).to(device)
