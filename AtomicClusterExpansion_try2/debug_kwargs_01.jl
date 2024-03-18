using MyACE1x

function debug_kwargs(; kwargs_in...)
    kwargs = MyACE1x._clean_args(kwargs_in)
    # kwargs is a NamedTuple
    for (k,v) in zip(keys(kwargs), kwargs)
        println(k, " => ", v)
    end
    return kwargs
end

kwargs = debug_kwargs(
    elements = [:Ti, :Al],
    order = 3,
    totaldegree = 6,
    rcut = 5.5,
    r0 = r0
)