Example processing for one configuration:
```pyconsole
>>> coords.shape
torch.Size([1, 22, 3])

>>> coords[:,None].shape
torch.Size([1, 1, 22, 3])

>>> coords[:,:,None].shape
torch.Size([1, 22, 1, 3])
```

Using `None` in indexing insert a "dummy" axis (singleton dimension)
similar to np.newaxis?
```pyconsole
>>> AA = np.random.rand(1, 4, 3);

>>> AA.shape
(1, 4, 3)

>>> AA[:,np.newaxis].shape
(1, 1, 4, 3)

>>> AA[:,:,np.newaxis].shape
(1, 4, 1, 3)
```



Computing pairwise distance

```pyconsole
>>> rij = coords[:, :, None] - coords[:, None]

>>> rij.shape
torch.Size([1, 22, 22, 3])

>>> rij[0,0,0,:]
tensor([0., 0., 0.])

>>> rij[0,0,1,:]
tensor([-1.0141,  0.3929, -0.0730])

>>> rij[0,1,0,:]
tensor([ 1.0141, -0.3929,  0.0730])

>>> rij[0,1,1,:]
tensor([0., 0., 0.])

>>> rij[0,1,2,:]
tensor([-0.2926,  0.4259,  0.9597])

>>> rij[0,2,1,:]
tensor([ 0.2926, -0.4259, -0.9597])
```



Compute the norm:
```pyconsole
>>> dij = torch.norm(rij, dim=3)
```


# After removing self-interaction

Example

```pyconsole
>>> rij_[0,1,1,:] # before removing self interaction
tensor([0., 0., 0.])

>>> rij[0,1,1,:] # after removing self interaction
Out[19]: tensor([-0.2926,  0.4259,  0.9597])

>>> rij_[0,1,2,:] # this is the value that get "shifted-up"
tensor([-0.2926,  0.4259,  0.9597])
```

