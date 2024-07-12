# Case of periodic, non-orthonal cell

using LinearAlgebra: norm, det, dot, inv, I, kron
using SparseArrays: sparse

include("create_D1_matrix.jl")
include("create_D2_matrix.jl")

function main()

BCx = :PERIODIC_BC

LatVecs = zeros(Float64, 3, 3) # !!! stored by rows !!!!
LatVecs[1,:] = [0.667124384994991, 0.741249316661101, 0.074124931666110]*6
LatVecs[2,:] = [0.074124931666110, 0.741249316661101, 0.667124384994991]*6
LatVecs[3,:] = [0.667124384994991, 0.074124931666110, 0.741249316661101]*6
# LatVecs scale: 6


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

# Set up transformation matrices for non orthogonal cells
if CELL_TYPE == :NON_ORTHOGONAL
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
end



Lx, Ly, Lz = 6.0, 6.0, 6.0
Nx, Ny, Nz = 15, 15, 15
FDn = 6 # finite difference order
# using FDn*2 + 1 points

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


end