# Case of periodic, non-orthonal cell

using LinearAlgebra: norm, det, dot, inv

include("create_D1_matrix.jl")
include("create_D2_matrix.jl")

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



Lx = 6.0
Nx = 15 # Number of discretization points
FDn = 6 # finite difference order
# using FDn*2 + 1 points

D2mat = create_D2_matrix_FD(FDn, Nx, :PERIODIC_BC)
D1mat = create_D1_matrix_FD(FDn, Nx, :PERIODIC_BC)

#=
using Plots: heatmap, theme
using PlotThemes
theme(:dark)
heatmap(D2mat, yflip=true, aspect_ratio=:equal)
=#


[DL11,DL22,DL33,DG1,DG2,DG3] = blochLaplacian_1d(S,[0 0 0]);
if S.cell_typ < 3
    S.Lap_std = S.lapc_T(1,1) * kron(speye(S.Nz),kron(speye(S.Ny),DL11))  +  S.lapc_T(2,2) * kron(speye(S.Nz),kron(DL22,speye(S.Nx))) + ...
                S.lapc_T(3,3) * kron(DL33,kron(speye(S.Ny),speye(S.Nx))) ;
    if (S.cell_typ == 2)
        MDL = S.lapc_T(1,2) * kron(speye(S.Nz),kron(DG2,DG1))  +  S.lapc_T(2,3) * kron(DG3,kron(DG2,speye(S.Nx))) + ...
              S.lapc_T(1,3) * kron(DG3,kron(speye(S.Ny),DG1)) ;
        S.Lap_std = S.Lap_std + MDL;
    end
