using LinearAlgebra: dot, norm, det
include("gen_lattice.jl")

# Here we prepare LatVecs

# We also need to decide whether we store the vectors in rows or columns
STORED_BY = :COLUMNS
@assert STORED_BY in [:ROWS, :COLUMNS]

a = 2.0
b = 3.0
c = 5.0
α = 90.0
β = 90.0
γ = 60.0
LatVecs = gen_lattice_triclinic( a, b, c, α, β, γ )
# Using PWDFT.jl convention LatVecs are stored by columns

# If
if STORED_BY == :ROWS
    LatVecs = Matrix(LatVecs')
end

# Transformation matrix
Tinv = zeros(Float64, 3, 3)
Tinv[1,:] = [ 1, cos(deg2rad(γ)), cos(deg2rad(α)) ]
Tinv[2,:] = [ cos(deg2rad(γ)), 1, cos(deg2rad(β)) ]
Tinv[3,:] = [ cos(deg2rad(α)), cos(deg2rad(β)), 1 ]
T = inv(Tinv)


# Check the cell typ (orthogonal or non-orthogonal)
SMALL = 1e-12
if STORED_BY == :ROWS
    cond12 = abs( dot(LatVecs[1,:], LatVecs[2,:]) ) > SMALL
    cond23 = abs( dot(LatVecs[2,:], LatVecs[3,:]) ) > SMALL
    cond31 = abs( dot(LatVecs[3,:], LatVecs[1,:]) ) > SMALL
else
    # :COLUMNS
    cond12 = abs( dot(LatVecs[:,1], LatVecs[:,2]) ) > SMALL
    cond23 = abs( dot(LatVecs[:,2], LatVecs[:,3]) ) > SMALL
    cond31 = abs( dot(LatVecs[:,3], LatVecs[:,1]) ) > SMALL
end


CELL_TYPE = :ORTHOGONAL
if any([cond12, cond23, cond31]) # or using || operator
   CELL_TYPE = :NON_ORTHOGONAL
end

println("CELL_TYPE = $(CELL_TYPE)")
println("STORED_BY = $(STORED_BY)")



# Normalize
UnitLatVecs = zeros(Float64, 3, 3)
if STORED_BY == :ROWS
    UnitLatVecs[1,:] = LatVecs[1,:]/norm(LatVecs[1,:])
    UnitLatVecs[2,:] = LatVecs[2,:]/norm(LatVecs[2,:])
    UnitLatVecs[3,:] = LatVecs[3,:]/norm(LatVecs[3,:])
else
    # :COLUMNS
    UnitLatVecs[:,1] = LatVecs[:,1]/norm(LatVecs[:,1])
    UnitLatVecs[:,2] = LatVecs[:,2]/norm(LatVecs[:,2])
    UnitLatVecs[:,3] = LatVecs[:,3]/norm(LatVecs[:,3])
end


# Jacobian
if STORED_BY == :ROWS
    JacobianUnitLatVecs = det(UnitLatVecs')
else
    # :COLUMNS
    JacobianUnitLatVecs = det(UnitLatVecs)
end
if JacobianUnitLatVecs < 0
    error("Volume is negative!")
end

# metric_T, Gradient and Laplacian transformation matrices
if STORED_BY == :ROWS
    metric_T = UnitLatVecs * UnitLatVecs'
    grad_T = inv(UnitLatVecs')
else
    # :COLUMNS
    metric_T = UnitLatVecs' * UnitLatVecs
    grad_T = inv(UnitLatVecs)
end

lapc_T = grad_T * grad_T'


# Will metric_T, lapc_T depends on whether LatVecs are stored by columns or row?


println("These two should be the same:")

println("lapc_T = ")
display(lapc_T)

println("T = ")
display(T)
