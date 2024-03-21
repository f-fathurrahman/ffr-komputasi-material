import MyACE1
import MyACE1x
import MyJuLIP

function create_kwargs(; kwargs_in...)
    # Clean kwargs and convert to NamedTuple
    return MyACE1x._clean_args(kwargs_in)
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


# Extract some information from kwargs
elements = kwargs[:elements]
cor_order = MyACE1x._get_order(kwargs)
Deg, maxdeg, maxn = MyACE1x._get_degrees(kwargs)

# Restrict to pure2b only
pure2b = kwargs[:pure2b]
@assert pure2b

rbasis = MyACE1x._radial_basis(kwargs)

println("Deg = ", Deg)
println("maxdeg = ", maxdeg)
println("cor_order = ", cor_order)
println("delete2b = ", kwargs[:delete2b])

rpibasis = MyACE1x.Pure2b.pure2b_basis(
    species = MyJuLIP.AtomicNumber.(elements),
    Rn=rbasis,
    D=Deg,
    maxdeg=maxdeg,
    order=cor_order,
    delete2b = kwargs[:delete2b]
);
