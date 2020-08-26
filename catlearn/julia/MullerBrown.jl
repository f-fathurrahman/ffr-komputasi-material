struct MullerBrown
    p1::Vector{Float64}
    p2::Vector{Float64}
    p3::Vector{Float64}
    p4::Vector{Float64}
end

function MullerBrown()
    p1 = [-200.0, -1.0,  0.0, -10.0,  1.0, 0.0]
    p2 = [-100.0, -1.0,  0.0, -10.0,  0.0, 0.5]
    p3 = [-170.0, -6.5, 11.0,  -6.5, -0.5, 1.5]
    p4 = [15.0, 0.7, 0.6, 0.7, -1.0, 1.0]
    return MullerBrown(p1, p2, p3, p4)
end

function calc_energy_forces!(mb::MullerBrown, atoms::Atoms)
    p1 = mb.p1
    p2 = mb.p2
    p3 = mb.p3
    p4 = mb.p4

    x = atoms.positions[1,1]
    y = atoms.positions[2,1]

    z1 = p1[1] * exp( (p1[2]*(x-p1[5])^2) + (p1[3]*(x-p1[5])*(y-p1[6])) + (p1[4]*(y-p1[6])^2) )
    z2 = p2[1] * exp( (p2[2]*(x-p2[5])^2) + (p2[3]*(x-p2[5])*(y-p2[6])) + (p2[4]*(y-p2[6])^2) )
    z3 = p3[1] * exp( (p3[2]*(x-p3[5])^2) + (p3[3]*(x-p3[5])*(y-p3[6])) + (p3[4]*(y-p3[6])^2) )
    z4 = p4[1] * exp( (p4[2]*(x-p4[5])^2) + (p4[3]*(x-p4[5])*(y-p4[6])) + (p4[4]*(y-p4[6])^2) )
    
    atoms.energy = (z1 + z2 + z3 + z4)/100.0 # Scale

    dx1 = z1 * (2*p1[2]*(x-p1[5])+p1[3]*(y-p1[6]))
    dx2 = z2 * (2*p2[2]*(x-p2[5])+p2[3]*(y-p2[6]))
    dx3 = z3 * (2*p3[2]*(x-p3[5])+p3[3]*(y-p3[6]))
    dx4 = z4 * (2*p4[2]*(x-p4[5])+p4[3]*(y-p4[6]))
    
    dy1 = z1 * (p1[3]*(x-p1[5])+2*p1[4]*(y-p1[6]))
    dy2 = z2 * (p2[3]*(x-p2[5])+2*p2[4]*(y-p2[6]))
    dy3 = z3 * (p3[3]*(x-p3[5])+2*p3[4]*(y-p3[6]))
    dy4 = z4 * (p4[3]*(x-p4[5])+2*p4[4]*(y-p4[6]))

    Fx = dx1 + dx2 + dx3 + dx4
    Fy = dy1 + dy2 + dy3 + dy4
    Fz = 0.0
    atoms.forces[1,1] = -Fx/100.0
    atoms.forces[2,1] = -Fy/100.0
    atoms.forces[3,1] = -Fz/100.0

    return
end

function calc_energy_forces(mb::MullerBrown, x, y)
    p1 = mb.p1
    p2 = mb.p2
    p3 = mb.p3
    p4 = mb.p4

    z1 = p1[1] * exp( (p1[2]*(x - p1[5])^2) + (p1[3]*(x - p1[5])*(y - p1[6])) + (p1[4]*(y - p1[6])^2) )
    z2 = p2[1] * exp( (p2[2]*(x - p2[5])^2) + (p2[3]*(x - p2[5])*(y - p2[6])) + (p2[4]*(y - p2[6])^2) )
    z3 = p3[1] * exp( (p3[2]*(x - p3[5])^2) + (p3[3]*(x - p3[5])*(y - p3[6])) + (p3[4]*(y - p3[6])^2) )
    z4 = p4[1] * exp( (p4[2]*(x - p4[5])^2) + (p4[3]*(x - p4[5])*(y - p4[6])) + (p4[4]*(y - p4[6])^2) )
    
    energy = (z1 + z2 + z3 + z4)/100.0 # Scale

    dx1 = z1 * ( 2*p1[2]*(x - p1[5]) + p1[3]*(y - p1[6]) )
    dx2 = z2 * ( 2*p2[2]*(x - p2[5]) + p2[3]*(y - p2[6]) )
    dx3 = z3 * ( 2*p3[2]*(x - p3[5]) + p3[3]*(y - p3[6]) )
    dx4 = z4 * ( 2*p4[2]*(x - p4[5]) + p4[3]*(y - p4[6]) )
    
    dy1 = z1 * ( p1[3]*(x - p1[5]) + 2*p1[4]*(y - p1[6]) )
    dy2 = z2 * ( p2[3]*(x - p2[5]) + 2*p2[4]*(y - p2[6]) )
    dy3 = z3 * ( p3[3]*(x - p3[5]) + 2*p3[4]*(y - p3[6]) )
    dy4 = z4 * ( p4[3]*(x - p4[5]) + 2*p4[4]*(y - p4[6]) )

    Fx = dx1 + dx2 + dx3 + dx4
    Fy = dy1 + dy2 + dy3 + dy4
    Fz = 0.0

    forces = zeros(Float64,3,1)
    forces[1,1] = -Fx/100.0
    forces[2,1] = -Fy/100.0
    forces[3,1] = -Fz/100.0

    return energy, forces
end