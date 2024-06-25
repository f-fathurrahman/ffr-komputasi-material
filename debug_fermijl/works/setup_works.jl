Pkg.activate("../")

using Revise

# Guard against multiple push!
!( "./MyFermi" in LOAD_PATH) && push!(LOAD_PATH, "./MyFermi")
!( "./MyMolecules" in LOAD_PATH) && push!(LOAD_PATH, "./MyMolecules")
!( "./MyGaussianBasis" in LOAD_PATH) && push!(LOAD_PATH, "./MyGaussianBasis")
