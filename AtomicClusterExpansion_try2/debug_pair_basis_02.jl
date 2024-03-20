# This is adapted from function `ACE1x._pair_basis` defined in defaults.jl

import MyACE1
import MyACE1x
import MyJuLIP

function create_kwargs(; kwargs_in...)
    # Clean kwargs and convert to NamedTuple
    return MyACE1x._clean_args(kwargs_in)
end

#
# Helper functions (originally sub-functions)
#
function _s2i(s, trans_pair)
    return MyJuLIP.z2i(trans_pair.zlist, MyJuLIP.AtomicNumber(s))
end

function _r_basis(s1, s2, penv, allr0, allrcut, maxn)
    # create envelope first
    _env = MyACE1.PolyEnvelope(penv, allr0[(s1, s2)], allrcut[(s1, s2)] )
    # then the polynomials
    return MyACE1.OrthPolys.transformed_jacobi_env(
        maxn, alltrans[(s1, s2)], _env, allrcut[(s1, s2)]
    )
end

function _x_basis(s1, s2, pin, pcut, maxn)
    # simply pass maxn, trans, rcut, pcut and pin
    return MyACE1.OrthPolys.transformed_jacobi(
        maxn, alltrans[(s1, s2)], allrcut[(s1, s2)];
        pcut = pcut, pin = pin
    )
end


#
# Main script here
#

kwargs = create_kwargs(
    elements = [:Ti, :Al],
    order = 3,
    totaldegree = 6,
    rcut = 5.5,
    r0 = 2.88
);

rbasis_type = kwargs[:pair_basis]
elements = kwargs[:elements]

@assert rbasis_type == :legendre
@info "\n Using legendre as radial basis (asserted)"

if kwargs[:pair_degree] == :totaldegree
    @info "\n Using totaldegres as pair_degree"
    @info "\n maxn will be determined from degrees"
    Deg, maxdeg, maxn = MyACE1x._get_degrees(kwargs)
#
elseif kwargs[:pair_degree] isa Integer
    @info "\n Using specified pair_degree"
    maxn = kwargs[:pair_degree]
else
    error("Cannot determine `maxn` for pair basis from information provided.")
end
@info "\n maxn = $(maxn)"
# maxn will be passed to _r_basis or _x_basis which will be passed
# to transformed_jacobi

# Determine all rcut
allrcut = MyACE1x._get_all_rcut(kwargs; _rcut = kwargs[:pair_rcut])
if allrcut isa Number
    allrcut = Dict([(s1, s2) => allrcut for s1 in elements, s2 in elements]...)
end
@info "\n allrcut = $(allrcut)"

# Determine all transforms
trans_pair = MyACE1x._transform(kwargs, transform = kwargs[:pair_transform])
#
alltrans = Dict([
    (s1, s2) => trans_pair.transforms[_s2i(s1, trans_pair), _s2i(s2, trans_pair)].t
    for s1 in elements, s2 in elements]...)
@info "\n alltrans = $(alltrans)"


# Determine all r0
allr0 = MyACE1x._get_all_r0(kwargs)
@info "\n allr0 = $(allr0)"


# Construct all rbases
envelope = kwargs[:pair_envelope]
@info "\n envolope = $(envelope)"
if envelope isa Tuple
    if envelope[1] == :x
        @info "\n Using x envelope"
        pin = envelope[2]
        pcut = envelope[3]
        rbases = [ _x_basis(s1, s2, pin, pcut, maxn) for s1 in elements, s2 in elements ];
    elseif envelope[1] == :r
        @info "\n Using r envelope"
        penv = envelope[2]
        rbases = [ _r_basis(s1, s2, penv, allr0, allrcut, maxn) for s1 in elements, s2 in elements ];
    end
end

pairB = MyACE1.PolyPairBasis(rbases, elements);
# pairB.J is essentially rbases
# only pairB.J is now an `SMatrix` instead of `Matrix`

@info "Script ended here"

# TODO: investigate pairB.J.data