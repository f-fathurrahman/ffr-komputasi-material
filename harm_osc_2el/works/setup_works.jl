Pkg.activate("../")

using Revise

# Guard against multiple push!
!( "./HarmOsc2El" in LOAD_PATH) && push!(LOAD_PATH, "./HarmOsc2El")
!( "./MyModule" in LOAD_PATH) && push!(LOAD_PATH, "./MyModule")