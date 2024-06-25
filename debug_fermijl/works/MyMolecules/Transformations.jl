"""
    MyMolecules.translate!(atoms::Molecule, r::Vector{T})

Translate atoms into position r.
"""
function translate!(atoms::Molecule, r::Vector{T}) where T
    for a in atoms
        for i = 1:3
            a.xyz[i] -= r[i]
        end
    end
end

"""
    MyMolecules.translate(atoms::Molecule)

Return a copy of atoms translated into position r.
"""
function translate(atoms::Molecule, r::Vector{T}) where T
    newatoms = deepcopy(atoms)
    translate!(newatoms, r)
    return newatoms
end

"""
    MyMolecules.Cn(v::Vector, n::Int)

Returns an n rotation matrix about vector v
"""
function Cn(v, n::Integer)
    θ = 2*π / n
    return rotation_matrix(v, θ)
end

"""
    MyMolecules.rotation_matrix(v, θ::AbstractFloat)

Returns a rotation matrix by θ about a rotation axis v
"""
function rotation_matrix(V, θ::AbstractFloat)
    cosθ = cos(θ)
    sinθ = sin(θ)
    v = normalize(V)
    a = [1,2,3]
    O = zeros(eltype(v), (3,3))
    O .+= 1 - cosθ
    for i = 1:3, j = 1:3
        if i == j
            O[i,i] *= v[i]^2
            O[i,i] += cosθ
        else
            O[i,j] *= v[i]*v[j]
            b = [i,j]
            C = [i for i in a if i ∉ b][1]
            if i < j
                O[i,j] += (-1)^(i+j) * v[C]*sinθ
            else
                O[i,j] += (-1)^(i+j-1) * v[C]*sinθ
            end
        end
    end
    return O
end

function rotation_matrixd(v, θ::AbstractFloat)
    r = deg2rad(θ)
    return rotation_matrix(v, r)
end

"""
    MyMolecules.reflection_matrix(v::Vector)

Returns a reflection matrix through the plane with normal vector v
"""
function reflection_matrix(V)
    O = zeros(eltype(V), (3,3))
    v = normalize(V)
    for i = 1:3, j = i:3
        if i == j
            O[i,i] = 1 - 2*v[i]^2
        else
            O[i,j] = -2 * v[i] * v[j]
            O[j,i] = O[i,j]
        end
    end
    return O
end

function σ(v)
    return reflection_matrix(v)
end

"""
    MyMolecules.Sn(v::Vector, n::int)

Returns Sn improper rotation about vector v
"""    
function Sn(v, n)
    return Cn(v, n) * σ(v)
end

"""
    MyMolecules.inversion_matrix()

Returns a diagonal matrix with -1 along the diagonal
"""
function inversion_matrix()
    a = zeros(Int8, (3,3))
    for i = 1:3
        a[i,i] = -1
    end    
    return a
end

function i()
    return inversion_matrix()
end
    
"""
    MyMolecules.transform!(A::Molecule, O::Array)

Transforms xyz coordinates of each atom in A by some operation O
"""
function transform!(atoms::Molecule, O::Array)
    for a in atoms
        a.xyz .= O * a.xyz
    end
end

"""
    MyMolecules.transform(A::Molecule, O::Array)

Transforms xyz coordinates of each atom in A by some operation O
Returns new list of atoms with transformed xyz coordinates
"""
function transform(atoms::Molecule, O::Array)
    newatoms = deepcopy(atoms)
    transform!(newatoms, O)
    return newatoms
end

"""
    MyMolecules.issame(A::Molecule, B::Molecule)

Checks to see if two arrays of atoms are equivalent under permutation
Returns true or false
"""
function isequivalent(A::Molecule, B::Molecule)
    h = []
    len = size(A,1)
    #println(A, "\n", B)
    for i = 1:len
        for j = 1:len
            if A[i].mass == B[j].mass
                zs = broadcast(abs, A[i].xyz - B[j].xyz)
                c =  isapprox(zs, [0.0,0.0,0.0], atol=tol)
                #println(zs, " : ", c)
                #if isapprox(broadcast(abs, A[i].xyz - B[j].xyz), [0.0,0.0,0.0], rtol=1E-5)
                if c
                    push!(h, i)
                    break
                end
            end
        end
    end
    return size(h,1) == len
end