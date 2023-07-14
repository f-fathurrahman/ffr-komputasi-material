import ACE1x
import ACE1
import JuLIP
using NamedTupleTools: namedtuple
using JuLIP: AtomicNumber

const MY_KW_DEFAULTS = Dict(
    :elements => nothing,
    :order => nothing,
    :totaldegree => nothing,
    :wL => 1.5,
    #
    :rin => 0.0,
    :r0 => :bondlen,
    :rcut => (:bondlen, 2.5),
    :transform => (:agnesi, 2, 4),
    :envelope => (:x, 2, 2),
    :rbasis => :legendre,
    #
    :pure2b => true,
    :delete2b => true,
    :pure => false,
    #
    :pair_rin => :rin,
    :pair_rcut => :rcut,
    :pair_degree => :totaldegree,
    :pair_transform => (:agnesi, 1, 3),
    :pair_basis => :legendre,
    :pair_envelope => (:r, 2),
    #
    :Eref => missing
)



const MY_KW_ALIASES = Dict(
    :N => :order,
    :species => :elements,
    :trans => :transform,
)




function my_clean_args(kwargs)
    dargs = Dict{Symbol, Any}()

    # Check for aliases
    for key in keys(kwargs)
        if haskey(MY_KW_ALIASES, key)
            dargs[MY_KW_ALIASES] = kwargs[key]
        else
            dargs[key] = kwargs[key]
        end
    end

    # Add defaults if kwargs does not contain them
    for key in keys(MY_KW_DEFAULTS)
        if !haskey(dargs, key)
            dargs[key] = MY_KW_DEFAULTS[key]
        end
    end

    # Make pair_rcut equal to rcut
    if dargs[:pair_rcut] == :rcut
        dargs[:pair_rcut] = dargs[:rcut]
    end

   return namedtuple(dargs)
end


function my_mb_ace_basis(kwargs)
    elements = kwargs[:elements]
    cor_order = ACE1x._get_order(kwargs)
    Deg, maxdeg, maxn = ACE1x._get_degrees(kwargs)
    rbasis = ACE1x._radial_basis(kwargs)
    pure2b = kwargs[:pure2b]

    if pure2b && kwargs[:pure]
        # error("Cannot use both `pure2b` and `pure` options.")
        @info("Option `pure = true` overrides `pure2b=true`")
        pure2b = false
    end

    if pure2b
        rpibasis = ACE1x.Pure2b.pure2b_basis(
            species = AtomicNumber.(elements),
            Rn=rbasis,
            D=Deg,
            maxdeg=maxdeg,
            order=cor_order,
            delete2b=kwargs[:delete2b])
    elseif kwargs[:pure]
        dirtybasis = ACE1.ace_basis(
            species = AtomicNumber.(elements),
            rbasis=rbasis,
            D=Deg,
            maxdeg=maxdeg,
            N = cor_order, )
        _rem = kwargs[:delete2b] ? 1 : 0
        rpibasis = ACE1x.Purify.pureRPIBasis(dirtybasis; remove = _rem)
    else
        rpibasis = ACE1.ace_basis(
            species = AtomicNumber.(elements),
            rbasis=rbasis,
            D=Deg,
            maxdeg=maxdeg,
            N = cor_order, )
    end

    return rpibasis
end



function my_pair_basis(kwargs)

    println()
    println("-------------------")
    println("ENTER my_pair_basis")
    println("-------------------")
    println()

    rbasis = kwargs[:pair_basis]
    elements = kwargs[:elements]

    println("rbasis = ", rbasis)
    println("elements = ", elements)

    if rbasis isa ACE1.ScalarBasis

        return rbasis

    elseif rbasis == :legendre

        if kwargs[:pair_degree] == :totaldegree
            Deg, maxdeg, maxn = ACE1x._get_degrees(kwargs)
        elseif kwargs[:pair_degree] isa Integer
            maxn = kwargs[:pair_degree]
        else
            error("Cannot determine `maxn` for pair basis from information provided.")
        end

        println("maxn = ", maxn)

        allrcut = ACE1x._get_all_rcut(kwargs; _rcut = kwargs[:pair_rcut])
        println("allrcut = ", allrcut)
        if allrcut isa Number
            allrcut = Dict([(s1, s2) => allrcut for s1 in elements, s2 in elements]...)
        end
        # allrcut is converted to dict
        println("allrcut = ", allrcut)

        trans_pair = ACE1x._transform(kwargs, transform = kwargs[:pair_transform])
        println("trans_pair = ", trans_pair)

        _s2i(s) = JuLIP.Potentials.z2i(trans_pair.zlist, AtomicNumber(s))
    
        alltrans = Dict([(s1, s2) => trans_pair.transforms[_s2i(s1), _s2i(s2)].t
                       for s1 in elements, s2 in elements]...)
        println("alltrans = ")
        for t in alltrans
            println(t)
        end

        allr0 = ACE1x._get_all_r0(kwargs)

        function _r_basis(s1, s2, penv)
            _env = ACE1.PolyEnvelope(penv, allr0[(s1, s2)], allrcut[(s1, s2)] )
            return ACE1x.transformed_jacobi_env(maxn, alltrans[(s1, s2)], _env, allrcut[(s1, s2)])
        end

        _x_basis(s1, s2, pin, pcut)  = ACE1x.transformed_jacobi(maxn, alltrans[(s1, s2)], allrcut[(s1, s2)];
                                             pcut = pcut, pin = pin)

        envelope = kwargs[:pair_envelope]
        if envelope isa Tuple
            if envelope[1] == :x
                pin = envelope[2]
                pcut = envelope[3]
                rbases = [ _x_basis(s1, s2, pin, pcut) for s1 in elements, s2 in elements ]
            elseif envelope[1] == :r
                penv = envelope[2]
                rbases = [ _r_basis(s1, s2, penv) for s1 in elements, s2 in elements ]
            end
        end
    end

    println()
    println("-------------------")
    println("EXIT my_pair_basis")
    println("-------------------")
    println()

    return ACE1.PolyPairBasis(rbases, elements)
end



function my_ace_basis(; kwargs...)
    
    println()
    println("----------------------")
    println("kwargs before cleaned:")
    println("----------------------")
    for k in keys(kwargs)
        println(k, " : ", kwargs[k])
    end
   
    kwargs = my_clean_args(kwargs)

    println()
    println("---------------------")
    println("kwargs after cleaned:")
    println("---------------------")
    for k in keys(kwargs)
        println(k, " : ", kwargs[k])
    end

    rpiB = my_mb_ace_basis(kwargs)
   
    pairB = my_pair_basis(kwargs)

    return JuLIP.MLIPs.IPSuperBasis([pairB, rpiB])
end

basis = my_ace_basis(
    elements = [:Si,],
    order = 3,   
    totaldegree = 10,
    rcut = 5.0
)
@show typeof(basis)