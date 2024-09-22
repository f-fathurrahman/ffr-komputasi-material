import torch

torch.manual_seed(1234)

Natoms = 4
Ndata = 1
coords = torch.rand(Ndata,Natoms,3)
print("Coordinates: ")
for ia in range(Natoms):
    print(f"{ia+1} {coords[0,ia,:]}")

# coords have shape (Ndata,Natoms,3)
# Ndata will be referred to num_batches
# Natoms is num_channels

num_batches, num_channels, _ = coords.size()
# num_channels: no. of atoms (?)

# computes pairwise distance vectors
rij_ = coords[:, :, None] - coords[:, None]
dij_ = torch.norm(rij_, dim=3)
# coords.shape = torch.Size([1, 4, 3])
# coords[:,:,None].shape = torch.Size([1, 4, 1, 3])
# coords[:,None].shape = torch.Size([1, 1, 4, 3])
# rij_.shape = torch.Size([1, 4, 4, 3])


mask = ~torch.eye(num_channels, dtype=torch.bool, device=coords.device) # remove self-interaction
# inverse of unit matrix of type bool
# all diagonal elements are False
# other elements are True

rij = torch.masked_select(rij_, mask.unsqueeze(2)).view(num_batches, num_channels, num_channels - 1, 3)
# rij contains the vectors last dimension is 3

dij = torch.masked_select(dij_, mask).view(num_batches, num_channels, num_channels - 1)
# dij is the magnitude of the vector

dij_inv = 1 / dij
dij2_inv = dij_inv * dij_inv

loc_env_r = dij_inv  # distances
loc_env_a = rij * dij2_inv.unsqueeze(3)  # vectors

print("Info")
for ia in range(Natoms):
    for ja in range(Natoms-1):
        print(f"{ia} {ja} dist={dij[0,ia,ja]:.6f} r={loc_env_r[0,ia,ja]:.6f} a={loc_env_a[0,ia,ja,:]}")
