Where are basic data structures defined?

Struct `S`: some fields are added and modified.

# Hybrid DFT
S.xc == 41 --> using hydbrid DFT (PBE0)

S.ACEFlag = 1

Use Fourier space for exx:

S.ExxMethod == 'FOURIER_SPACE'
S.exxmethod = 0


# Laplacian setup

In setup_defaults.m

```matlab
S = lapIndicesValues_1d(S);
S = gradIndicesValues(S);

% Calculate discrete laplacian
[DL11,DL22,DL33,DG1,DG2,DG3] = blochLaplacian_1d(S,[0 0 0]);
% This is for S.cell_typ < 3
S.Lap_std = S.lapc_T(1,1) * kron(speye(S.Nz),kron(speye(S.Ny),DL11)) + ...
            S.lapc_T(2,2) * kron(speye(S.Nz),kron(DL22,speye(S.Nx))) + ...
            S.lapc_T(3,3) * kron(DL33,kron(speye(S.Ny),speye(S.Nx)));
if( S.cell_typ == 2 )
    MDL = S.lapc_T(1,2) * kron(speye(S.Nz),kron(DG2,DG1)) + ...
          S.lapc_T(2,3) * kron(DG3,kron(DG2,speye(S.Nx))) + ...
          S.lapc_T(1,3) * kron(DG3,kron(speye(S.Ny),DG1));
    S.Lap_std = S.Lap_std + MDL;
end
```




# Struct S

Some important fields of `S`:

- S.lat_vec

- S.cell_typ

- S.metric_T

- S.dx, S.dy, S.dz

- S.Nx, S.Ny, S.Nz

- S.Lap_std: Laplacian matrix?

- S.Atm(1), S.Atm(2), ...

