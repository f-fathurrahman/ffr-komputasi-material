mutable struct Electrons
    Nelectrons::Int64
    Nstates::Int64
    Nstates_occ::Int64
    Focc::Array{Float64,1}
    energies::Array{Float64,1}
end

function Electrons( Nelectrons::Int64; Nstates_extra=0 )

    is_odd = (Nelectrons%2 == 1)

    Nstates_occ = round(Int64, Nelectrons/2)
    if is_odd
        Nstates_occ = Nstates_occ + 1
    end
    
    Nstates = Nstates_occ + Nstates_extra
    
    Focc = zeros(Float64,Nstates)
    energies = zeros(Float64,Nstates)

    if !is_odd
        for i in 1:Nstates_occ
            i
        end
    else
    for i in 1:Nstates_occ
end