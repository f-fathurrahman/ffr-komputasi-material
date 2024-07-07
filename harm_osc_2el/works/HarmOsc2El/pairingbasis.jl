struct Pairing_State
    p::Int64
    s_z::Int64
    g::Float64
end

function Pairing_State(label, g)
    if label%2 != 0
        s_z = 1
        p = (label-1) / 2 + 1
    else
        s_z = -1
        p = (label-2) / 2 + 1
    end

    return Pairing_State(p, s_z, g)
end

struct Pairing <: Basis
    l::Int64
    states::Vector{Pairing_State}
end

function Pairing(l::Int64, g::Float64)
    states = [Pairing_State(label, g) for label in 1:l]
    
    return Pairing(l, states)
end

function Ĥ₀(p::T, q::T) where T <: Pairing_State
    """
    The one body part of the hamiltonian
    
    delta = 1
    """
    if p == q # if p and q are the same state
        return p.p - 1
    else
        return 0.0
    end
end

function V̂(p::T, q::T, r::T, s::T) where T <: Pairing_State
    g = p.g
    p1, s1 = p.p, p.s_z
    p2, s2 = q.p, q.s_z
    p3, s3 = r.p, r.s_z
    p4, s4 = s.p, s.s_z

    if p1 != p2 || p3 != p4
        return 0.0
    end
    if s1 == s2 || s3 == s4
        return 0.0
    end
    if s1 == s3 && s2 == s4
        return -g/2
    end
    if s1 == s4 && s2 == s3
        return g/2
    end
end

function pairing_exact(g)
    H = [ 2-g  -g/2 -g/2 -g/2 -g/2  0
         -g/2  4-g  -g/2 -g/2  0   -g/2
         -g/2 -g/2   6-g  0   -g/2 -g/2
         -g/2 -g/2   0    6-g -g/2 -g/2
         -g/2  0    -g/2 -g/2  8-g -g/2
          0   -g/2  -g/2 -g/2 -g/2  10-g]
    return la.eigvals(H)[1]
end

function pairing_MBPT2(g)
    return -g^2 / 4 * ( 1 / (4+g) + 1/(6+g) + 1/(2+g) + 1/(4+g) ) + (2 - g)
end
;