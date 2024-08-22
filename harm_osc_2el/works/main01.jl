using HarmOsc2El

function setup_rhf_sys(;n, l)
    ω = 0.25
    basis = SpinBasis(HOBasis(l, ω))

    V = HOCoulomb(ω, shielding = 0.25)
    grid = [x for x in range(-10, stop = 10, length = 2001)]
    system = System(n, basis, grid, V);

    rhf = RHF(system)
    t = @elapsed compute_ground_state!(rhf);
    rhf_system = System(rhf)
    println("Reference energy: $(reference_energy(rhf_system)), $(t) s")

    return rhf_system
end

# rename?
function hf_energy(; n, l)
    ω = 0.25
    basis = SpinBasis(HOBasis(l, ω))

    V = HOCoulomb(ω, shielding = 0.25)
    grid = [x for x in range(-10, stop = 10, length = 2001)]
    system = System(n, basis, grid, V);

    hf = HF(system)
    t = @elapsed compute_ground_state!(hf);
    println("HF energy: $(energy(hf)), $(t) s")
end

hf_energy(; n=2, l=20)
