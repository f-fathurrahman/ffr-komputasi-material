import MyACE1
import MyACE1x
import MyJuLIP

function debug_pair_basis(; kwargs_in...)

    # Clean kwargs and convert to NamedTuple
    kwargs = MyACE1x._clean_args(kwargs_in)


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
    _s2i(s) = MyJuLIP.z2i(trans_pair.zlist, MyJuLIP.AtomicNumber(s))
    alltrans = Dict([(s1, s2) => trans_pair.transforms[_s2i(s1), _s2i(s2)].t
                     for s1 in elements, s2 in elements]...)

    # Determine all r0
    allr0 = MyACE1x._get_all_r0(kwargs)

    function _r_basis(s1, s2, penv)
        _env = MyACE1.PolyEnvelope(penv, allr0[(s1, s2)], allrcut[(s1, s2)] )
        return MyACE1.OrthPolys.transformed_jacobi_env(
            maxn, alltrans[(s1, s2)], _env, allrcut[(s1, s2)]
        )
    end

    function _x_basis(s1, s2, pin, pcut)
        return MyACE1.OrthPolys.transformed_jacobi(
            maxn, alltrans[(s1, s2)], allrcut[(s1, s2)];
            pcut = pcut, pin = pin
        )
    end

    # Construct all rbases
    envelope = kwargs[:pair_envelope]
    if envelope isa Tuple
        if envelope[1] == :x
            @info "Pass here 60"
            pin = envelope[2]
            pcut = envelope[3]
            rbases = [ _x_basis(s1, s2, pin, pcut) for s1 in elements, s2 in elements ];
        elseif envelope[1] == :r
            @info "Pass here 65"
            penv = envelope[2]
            rbases = [ _r_basis(s1, s2, penv) for s1 in elements, s2 in elements ];
        end
    end
    
    @info "Pass here 71"

    println("typeof rbases = ", typeof(rbases))
    println("typeof elements = ", typeof(elements))

    return MyACE1.PolyPairBasis(rbases, elements)
end


pairB = debug_pair_basis(
    elements = [:Ti, :Al],
    order = 3,
    totaldegree = 6,
    rcut = 5.5,
    r0 = 2.88
);

@info "Pass here 92" # to prevent print out?

