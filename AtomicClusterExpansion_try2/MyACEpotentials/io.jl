using Pkg

export load_potential
export save_potential

"""
    function load_potential(fname::AbstractString; new_format=false, verbose=true)

Load MyACE potential from given file `fname`.

# Kwargs
- `new_format=false` - If true returns potential as `MyACEmd.MyACEpotential` format, else use old MyJuLIP format
- `verbose=true`     - Display version info on load
"""
function load_potential(
    fname::AbstractString;
    new_format=false,
    verbose=true
)
    pot_tmp = load_dict(fname)
    if verbose && haskey(pot_tmp, "Versions")
        println("\nThis potential was saved with following versions:\n")
        for (k,v) in pot_tmp["Versions"]
            n = VersionNumber(v["major"], v["minor"], v["patch"])
            println(k," v",n)
        end
        println("\n", "If you have problems with using this potential, pin your installation to above versions.\n")
    end
    if haskey(pot_tmp, "IP")
        pot = read_dict(pot_tmp["IP"])
    elseif haskey(pot_tmp, "potential")
        pot = read_dict(pot_tmp["potential"])
    else
        error("Potential format not recognised")
    end
    if new_format
        return MyACEpotential(pot.components)
    else
        return pot
    end
end


"""
    save_potential( fname, potential::MyACE1x.MyACE1Model; save_version_numbers=true, meta=nothing)

Save MyACE potentials. Prefix is either .json, .yml or .yace, which also determines file format.

# Kwargs
- save_version_numbers=true  : If true save version information or relevant packages
- `meta=nothing`             : Seve some metadata with the potential (needs to be `Dict{String, Any}`)
"""
function save_potential( fname, potential::MyACE1x.MyACE1Model; save_version_numbers=true, meta=nothing)
    return save_potential(fname, potential.potential; save_version_numbers=save_version_numbers, meta=meta)
end

function save_potential( fname, potential::MyACEmd.MyACEpotential; save_version_numbers=true, meta=nothing)
    return save_potential(fname, potential.potentials; save_version_numbers=save_version_numbers, meta=meta)
end

function save_potential(fname, potential; save_version_numbers=true, meta=nothing)
    if save_version_numbers
        versions = Dict()
        versions["MyACEpotentials"] = extract_version("MyACEpotentials")
        versions["MyACE1"] = extract_version("MyACE1")
        versions["MyACE1x"] = extract_version("MyACE1x")
        versions["MyACEmd"] = extract_version("MyACEmd")
        versions["MyACEbase"] = extract_version("MyACEbase")
        versions["MyJuLIP"] = extract_version("MyJuLIP")
        versions["MyACEfit"] = extract_version("MyACEfit")

        data = Dict(
            "IP" => write_dict(potential),
            "Versions" => versions
        )
    else
        data = Dict(
            "IP" => write_dict(potential)
        )
    end
    if !isnothing(meta)
        @assert isa(meta, Dict{String, <:Any}) "meta needs to be a Dict{String, Any}"
        data["meta"] = convert(Dict{String, Any}, meta)
    end
    save_dict(fname, data)
end


# used to extraction version numbers when saving
function extract_version(name::AbstractString)
    vals = Pkg.dependencies()|> values |> collect
    hit = filter(x->x.name==name, vals) |> only
    return hit.version
end


## Deprecations

@deprecate export2json(fname, model; meta=nothing) save_potential(fname, model; meta=meta)

