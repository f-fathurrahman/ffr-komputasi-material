using StaticArrays


function redirect_to_files(dofunc, outfile, errfile)
    open(outfile, "w") do out
        open(errfile, "w") do err
            redirect_stdout(out) do
                redirect_stderr(err) do
                    dofunc()
                end
            end
        end
    end
end

function show_fields(a)
    t = typeof(a)
    println("<div> ENTER Type of variable: ", t)
    for s in fieldnames(t)
        f = getfield(a, s)
        println("\n  fieldname = $s, type = $(typeof(f))")
    end
    println("</div> EXIT Type of variable: ", t)
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

function _is_collection_type(t)
    is_Tuple = t <: Tuple
    is_SMatrix = t <: SMatrix
    is_SVector = t <: SVector
    is_Vector = t <: Vector
    is_Matrix = t <: Matrix
    return is_Tuple || is_SMatrix || is_SVector || is_Vector || is_Matrix
end


function show_fields_recursive(a_in; indent="  ")
    if _is_collection_type(typeof(a_in))
        #println("\nThis is a collection type: using its first element instead\n")
        a = first(a_in)
    else
        a = a_in
    end
    t = typeof(a)
    if t <: Number
        #@info "Stopped here because t=$t is a Number"
        #println(indent * "EXIT type: $t")
        return
    end
    if t <: UnitRange
        #println(indent * "EXIT type: $t")
        return
    end
    if t <: MyJuLIP.Potentials.SZList
        #println(indent * "EXIT type: $t")
        return
    end

    println()
    println(indent * "<div> ENTER type: $t")
    
    if t <: Dict
        #println(indent * "!!! This is a Dict: we will skip it")
        #println(indent * "EXIT type: $t")
        return
    end


    for s in fieldnames(t)
        f = getfield(a, s)
        println("\n" * indent * " - fieldname = $s, type = $(typeof(f))")
        show_fields_recursive(f, indent=indent*"  ")
    end

    println(indent * "</div> EXIT type: $t")
    return
end
