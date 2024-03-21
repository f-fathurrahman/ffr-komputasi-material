function show_fields(a)
    t = typeof(a)
    println("\n Type of variable: ", t)
    for s in fieldnames(t)
        f = getfield(a, s)
        println("\n  fieldname = $s, type = $(typeof(f))")
    end
end