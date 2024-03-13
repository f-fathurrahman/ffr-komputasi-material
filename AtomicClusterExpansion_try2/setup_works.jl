Pkg.activate("ACESUITE", shared=true)

using Revise

# Guard against multiple push!
!( "./MyJuLIP" in LOAD_PATH) && push!(LOAD_PATH, "./MyJuLIP")
!("./MyACEbase" in LOAD_PATH) && push!(LOAD_PATH, "./MyACEbase")
!("./MyACE1" in LOAD_PATH) && push!(LOAD_PATH, "./MyACE1")
!("./MyACE1x" in LOAD_PATH) && push!(LOAD_PATH, "./MyACE1x")
!("./MyACEfit" in LOAD_PATH) && push!(LOAD_PATH, "./MyACEfit")
!("./MyACEmd" in LOAD_PATH) && push!(LOAD_PATH, "./MyACEmd")
!("./MyACEpotentials" in LOAD_PATH) && push!(LOAD_PATH, "./MyACEpotentials")
