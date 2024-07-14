using LinearAlgebra: dot, norm, det
include("gen_lattice.jl")

BCx = :PERIODIC_BC
STORED_BY = :COLUMNS

a = 2.0
b = 3.0
c = 5.0
α = 90.0
β = 90.0
γ = 60.0

# Using PWDFT.jl convention LatVecs are stored by columns
#LatVecs = gen_lattice_fcc(6.0)
LatVecs = gen_lattice_triclinic( a, b, c, α, β, γ )

if STORED_BY == :COLUMNS
    LatVecs = LatVecs'
    # transpose because the code below assumes stored by rows
end


# Transformation matrix
Tinv = zeros(Float64, 3, 3)
Tinv[1,:] = [ 1, cos(deg2rad(γ)), cos(deg2rad(α)) ]
Tinv[2,:] = [ cos(deg2rad(γ)), 1, cos(deg2rad(β)) ]
Tinv[3,:] = [ cos(deg2rad(α)), cos(deg2rad(β)), 1 ]
T = inv(Tinv)


# Check the cell typ (orthogonal or non-orthogonal)
SMALL = 1e-12
cond12 = abs( dot(LatVecs[1,:], LatVecs[2,:]) ) > SMALL
cond23 = abs( dot(LatVecs[2,:], LatVecs[3,:]) ) > SMALL
cond31 = abs( dot(LatVecs[3,:], LatVecs[1,:]) ) > SMALL
CELL_TYPE = :ORTHOGONAL
if any([cond12, cond23, cond31]) # or using || operator
   CELL_TYPE = :NON_ORTHOGONAL
end
println("CELL_TYPE = $(CELL_TYPE)")



# Normalize
UnitLatVecs = zeros(Float64, 3, 3)
UnitLatVecs[1,:] = LatVecs[1,:]/norm(LatVecs[1,:])
UnitLatVecs[2,:] = LatVecs[2,:]/norm(LatVecs[2,:])
UnitLatVecs[3,:] = LatVecs[3,:]/norm(LatVecs[3,:])

# Jacobian
JacobianUnitLatVecs = det(UnitLatVecs')
if JacobianUnitLatVecs < 0
    error("Volume is negative!")
end
# metric_T, Gradient and Laplacian transformation matrices
metric_T = UnitLatVecs * UnitLatVecs'
#metric_T[1,2] = 2*metric_T[1,2]  # factor of 2 removed?
#metric_T[2,3] = 2*metric_T[2,3] 
#metric_T[1,3] = 2*metric_T[1,3]
grad_T = inv(UnitLatVecs')
lapc_T = grad_T * grad_T'
#lapc_T[1,2] = 2*lapc_T[1,2]
#lapc_T[2,3] = 2*lapc_T[2,3]
#lapc_T[1,3] = 2*lapc_T[1,3]

# Will metric_T, lapc_T depends on whether LatVecs are stored by columns or row?


println("lapc_T")
display(lapc_T); println()
# These are used for Laplacian matrix
println("lapc_T[1,2] = ", lapc_T[1,2])
println("lapc_T[2,3] = ", lapc_T[2,3])
println("lapc_T[1,3] = ", lapc_T[1,3])
