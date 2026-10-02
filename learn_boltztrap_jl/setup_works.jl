Pkg.activate("BOLTZTRAP", shared=true)

using Revise

# Guard against multiple push!
!( "./MyBoltzTraP/src" in LOAD_PATH) && push!(LOAD_PATH, "./MyBoltzTraP/src")
