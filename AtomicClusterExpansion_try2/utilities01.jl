function show_fields(a)
    t = typeof(a)
    println("\n Type of variable: ", t)
    for s in fieldnames(t)
        f = getfield(a, s)
        println("\n  fieldname = $s, type = $(typeof(f))")
    end
    return
end

function _is_type_skipped(t; indent=" ")
    if t <: Tuple
        println(indent * "Stopped here")
        return true
    end
    if t <: StaticArrays.SMatrix
        println(indent * "Stopped here")
        return true
    end
    if t <: StaticArrays.SVector
        println(indent * "Stopped here")
        return true
    end
    if t <: Vector
        println(indent * "Stopped here")
        return true
    end
    return false
end

function show_fields_recursive(a; indent=" ")
    t = typeof(a)
    println()
    println(indent * "------------------------")
    println(indent * "Parent type: $t")
    println(indent * "------------------------")
    
    if _is_type_skipped(t, indent=indent)
        return
    end

    for s in fieldnames(t)
        f = getfield(a, s)
        println("\n" * indent * " - fieldname = $s, type = $(typeof(f))")
    end
    #
    for s in fieldnames(t)
        f = getfield(a, s)
        t = typeof(f)
        if _is_type_skipped(t, indent=indent)
            continue
        end
        if isbitstype(t)
            continue
        end
        if t <: Vector
            @info "a is a Vector, processing its element instead"
            a = first(a)
        end
        show_fields_recursive(f, indent=indent*" ")
    end
    return
end
