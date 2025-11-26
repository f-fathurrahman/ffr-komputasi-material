Pkg.activate("../")

using Revise

# Guard against multiple push!
!( "./MyDFTK" in LOAD_PATH) && push!(LOAD_PATH, "./MyDFTK")
