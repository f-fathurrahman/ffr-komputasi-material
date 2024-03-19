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
    _env = MyACE1.PolyEnvelope(penv, allr0[(s1, s2)], allrcut[(s1, s2)] )
    return MyACE1.OrthPolys.transformed_jacobi_env(
        maxn, alltrans[(s1, s2)], _env, allrcut[(s1, s2)]
    )
end

function _x_basis(s1, s2, pin, pcut, maxn)
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

rbasis = kwargs[:pair_basis]
elements = kwargs[:elements]

@assert rbasis == :legendre
@info "Using legendre as radial basis (asserted)"

if kwargs[:pair_degree] == :totaldegree
    println("Pass here 19")
    Deg, maxdeg, maxn = MyACE1x._get_degrees(kwargs)
#
elseif kwargs[:pair_degree] isa Integer
    maxn = kwargs[:pair_degree]
#
else
    error("Cannot determine `maxn` for pair basis from information provided.")
end
@info "maxn = $maxn"

# Determine all rcut
allrcut = MyACE1x._get_all_rcut(kwargs; _rcut = kwargs[:pair_rcut])
if allrcut isa Number
    allrcut = Dict([(s1, s2) => allrcut for s1 in elements, s2 in elements]...)
end

# Determine all transforms
trans_pair = MyACE1x._transform(kwargs, transform = kwargs[:pair_transform])
#
alltrans = Dict([
    (s1, s2) => trans_pair.transforms[_s2i(s1, trans_pair), _s2i(s2, trans_pair)].t
    for s1 in elements, s2 in elements]...)

# Determine all r0
allr0 = MyACE1x._get_all_r0(kwargs)

# Construct all rbases
envelope = kwargs[:pair_envelope]
if envelope isa Tuple
    if envelope[1] == :x
        @info "Pass here 60"
        pin = envelope[2]
        pcut = envelope[3]
        rbases = [ _x_basis(s1, s2, pin, pcut, maxn) for s1 in elements, s2 in elements ];
    elseif envelope[1] == :r
        @info "Pass here 65"
        penv = envelope[2]
        rbases = [ _r_basis(s1, s2, penv, allr0, allrcut, maxn) for s1 in elements, s2 in elements ];
    end
end

pairB = MyACE1.PolyPairBasis(rbases, elements);

@info "Script ended here"

# TODO: investigate pairB.J.data