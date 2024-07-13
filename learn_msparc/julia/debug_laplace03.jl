# Case of periodic, non-orthonal cell

using LinearAlgebra: norm, det, dot, inv, I, kron
using SparseArrays: sparse

include("create_D1_matrix.jl")
include("create_D2_matrix.jl")

BCx = :PERIODIC_BC

LatVecs = zeros(Float64, 3, 3) # !!! stored by rows !!!!
LatVecs[1,:] = [10.0, 0.0, 0.0]
LatVecs[2,:] = [0.0, 12.0, 0.0]
LatVecs[3,:] = [0.0, 1.0, 10.0]

Lx = norm(LatVecs[1,:])
Ly = norm(LatVecs[2,:])
Lz = norm(LatVecs[3,:])
Nx = 26
Ny = 28
Nz = 26
FDn = 6 # finite difference order
# using FDn*2 + 1 points


# Check the cell typ (orthogonal or non-orthogonal)
SMALL = 1e-12
cond12 = abs( dot(LatVecs[1,:], LatVecs[2,:]) ) > SMALL
cond23 = abs( dot(LatVecs[2,:], LatVecs[3,:]) ) > SMALL
cond31 = abs( dot(LatVecs[3,:], LatVecs[1,:]) ) > SMALL
CELL_TYPE = :ORTHOGONAL
if any([cond12, cond23, cond31]) # or using || operator
   CELL_TYPE = :NON_ORTHOGONAL
end

# Normalize
UnitLatVecs = zeros(Float64, 3, 3)
UnitLatVecs[1,:] = LatVecs[1,:]/norm(LatVecs[1,:])
UnitLatVecs[2,:] = LatVecs[2,:]/norm(LatVecs[2,:])
UnitLatVecs[3,:] = LatVecs[3,:]/norm(LatVecs[3,:])

# Set up transformation matrices

# Jacobian
JacobianUnitLatVecs = det(UnitLatVecs')
if JacobianUnitLatVecs < 0
    error("Volume is negative!")
end
# metric_T, Gradient and Laplacian transformation matrices
metric_T = UnitLatVecs * UnitLatVecs'
metric_T[1,2] = 2*metric_T[1,2] 
metric_T[2,3] = 2*metric_T[2,3] 
metric_T[1,3] = 2*metric_T[1,3]
grad_T = inv(UnitLatVecs')
lapc_T = grad_T * grad_T'
lapc_T[1,2] = 2*lapc_T[1,2]
lapc_T[2,3] = 2*lapc_T[2,3]
lapc_T[1,3] = 2*lapc_T[1,3]


DL11 = create_D2_matrix_FD(FDn, Lx, Nx, :PERIODIC_BC)
DG1 = create_D1_matrix_FD(FDn, Lx, Nx, :PERIODIC_BC)

DL22 = create_D2_matrix_FD(FDn, Ly, Ny, :PERIODIC_BC)
DG2 = create_D1_matrix_FD(FDn, Ly, Ny, :PERIODIC_BC)

DL33 = create_D2_matrix_FD(FDn, Lz, Nz, :PERIODIC_BC)
DG3 = create_D1_matrix_FD(FDn, Lz, Nz, :PERIODIC_BC)

#=
using Plots: heatmap, theme
using PlotThemes
theme(:dark)
heatmap(D2mat, yflip=true, aspect_ratio=:equal)
=#

speye(N::Int64) = sparse(Matrix(1.0I, N, N))

Ix = speye(Nx)
Iy = speye(Ny)
Iz = speye(Nz)

Lap_std = lapc_T[1,1] * kron( Iz, kron(Iy, DL11) ) + 
          lapc_T[2,2] * kron( Iz, kron(DL22, Ix) ) +
          lapc_T[3,3] * kron( DL33, kron(Iy, Ix) )
if CELL_TYPE == :NON_ORTHOGONAL
    MDL = lapc_T[1,2] * kron( Iz, kron(DG2, DG1) ) +
          lapc_T[2,3] * kron( DG3, kron(DG2, Ix) ) +
          lapc_T[1,3] * kron( DG3, kron(Iy, DG1) )
    Lap_std += MDL
end

Npoints = Nx*Ny*Nz
Nstates = 3

X = zeros(Float64, Npoints, Nstates)
X[100,1] = 1.0
X[100,2] = 2.0
X[100,3] = 3.0

# Apply Laplacian to X

if CELL_TYPE == :ORTHOGONAL
    X1 = reshape(X, Nx, Ny, :)
    l_zs = size(X1, 3)
    Hlx1 = zeros(Float64, Nx, Ny, l_zs)
    for i in 1:l_zs
        @views Hlx1[:,:,i] = DL11*X1[:,:,i] + X1[:,:,i]*DL22'
    end
    Hlx1 = reshape(Hlx1, Nx*Ny, Nz, :)
    X2 = reshape(X, Nx*Ny, Nz, :)
    l_s = size(X2, 3)
    Hlx2 = zeros(Nx*Ny, Nz, l_s)
    for i in 1:l_s
        @views Hlx2[:,:,i] = X2[:,:,i]*DL33'
    end    
    DLX = reshape(Hlx1 + Hlx2, :, l_s)
else
    T11 = lapc_T[1,1]
    T22 = lapc_T[2,2]
    T33 = lapc_T[3,3]
    T12 = lapc_T[1,2]
    T23 = lapc_T[2,3]
    T13 = lapc_T[1,3]
    #
    X1 = reshape(X, Nx, Ny, :)
    l_zs = size(X1, 3)
    Hlx1 = zeros(Nx, Ny, l_zs)
    Hlx2 = zeros(Nx, Ny, l_zs)
    for i in 1:l_zs
        @views Hlx1[:,:,i] = T11*DL11*X1[:,:,i] + T22*X1[:,:,i]*DL22' + T12*DG1*X1[:,:,i]*DG2'
        @views Hlx2[:,:,i] = T13*DG1*X1[:,:,i] + T23*X1[:,:,i]*DG2'
    end
    Hlx1 = reshape(Hlx1, Nx*Ny, Nz, :)
    Hlx2 = reshape(Hlx2, Nx*Ny, Nz, :)
    X2 = reshape(X, Nx*Ny, Nz, :)
    l_s = size(X2, 3)
    for i in 1:l_s
        @views Hlx2[:,:,i] = T33*X2[:,:,i]*DL33' + Hlx2[:,:,i]*DG3'
    end
    DLX = reshape(Hlx1 + Hlx2, :, l_s)
end


