using LinearAlgebra: dot, norm
using Printf
using SparseArrays
using SpecialFunctions: besselj
using Infiltrator

function gen_idx_xy2ip(Nx, Ny)
    idx_xy2ip = zeros(Int64, Nx, Ny)
    ip = 0
    for i in 1:Nx, j in 1:Ny # the same ordering as in Python
        ip += 1
        idx_xy2ip[i,j] = ip
    end
    return idx_xy2ip
end


function init_H_and_V(Nx, Ny, W)

    idx_xy2ip = gen_idx_xy2ip(Nx, Ny)

    Npoints = Nx * Ny # total number of sites

    row_H = zeros(Int64, Npoints*4 - Nx*2) # the row indices for H
    col_H = zeros(Int64, Npoints*4 - Nx*2) # the column indices for H
    Hij = -ones(ComplexF64, Npoints*4 - Nx*2) # nonzero Hamiltonian elements

    row_V = zeros(Int64, Npoints*2) # row indices for V
    col_V = zeros(Int64, Npoints*2) # column indices for V
    Vij = zeros(ComplexF64, Npoints*2) # nonzero velocity matrix elements

    row_U = 1:Npoints # row and column indices for U
    col_U = 1:Npoints
    Uij = rand(Npoints)*W .- W/2

    count_H = 0
    count_V = 0

    for i in 1:Nx, j in 1:Ny # the same ordering as in Python
        #
        # (0) # get the index of the center site
        idx_center = idx_xy2ip[i,j]
        #
        # (1) consider the left neighbor (periodic boundary)
        ii = i - 1
        if ii <= 0
            ii += Nx # wrap
        end
        idx_left = idx_xy2ip[ii,j]
        #
        count_H += 1    
        row_H[count_H] = idx_center
        col_H[count_H] = idx_left
        #
        count_V += 1
        row_V[count_V] = idx_center
        col_V[count_V] = idx_left
        Vij[count_V] = im
        #
        # (2) consider the right neighbor (periodic boundary)
        ii = i + 1
        if ii >= Nx + 1
            ii -= Nx # wrap
        end
        idx_right = idx_xy2ip[ii,j]
        #
        count_H += 1
        row_H[count_H] = idx_center
        col_H[count_H] = idx_right
        #
        count_V += 1
        row_V[count_V] = idx_center
        col_V[count_V] = idx_right
        Vij[count_V] = -im
        #
        # (3) consider the upper neighbor (open boundary)
        #@printf("\n -- i=%3d j=%3d idx_center=%d\n", i, j, idx_center)
        #@printf("Periodic BC: idx_left=%3d idx_right=%3d\n" , idx_left, idx_right)
        if j < Ny
            idx_up = idx_xy2ip[i,j+1]
            count_H += 1
            row_H[count_H] = idx_center
            col_H[count_H] = idx_up
            #@printf("Open BC: idx_up=%3d\n", idx_up)
        end
        #
        # (4) consider the down neighbor (open boundary)
        if j > 1
            idx_down = idx_xy2ip[i,j-1]
            count_H += 1
            row_H[count_H] = idx_center
            col_H[count_H] = idx_down
            #@printf("Open BC: index_down=%3d\n", idx_down)
        end
    end

    H = sparse( row_H, col_H, Hij )
    U = sparse( row_U, col_U, Uij )
    H = H + U
    V = sparse( row_V, col_V, Vij )

    return H, V

end # function

function init_state(N)
    θ = 2π * rand(N)
    ϕ = exp.(im*θ)
    ϕ ./= norm(ϕ);
    return ϕ
end

function jackson_damping(M)
    m = 1:M
    a = 1.0 / (M + 1.0)
    g = (1.0 .- m * a) .* cos.(pi * m * a)
    g .+= a * sin.(pi * m * a) / tan.(pi * a)
    return g
end

function calc_moments(M, H_scaled, ϕ_l, ϕ_r)
    #
    C = zeros(Float64, M)
    #
    ϕ_0 = copy(ϕ_r)
    C[1] = real(dot(ϕ_l, ϕ_0))
    #
    ϕ_1 = H_scaled * ϕ_0
    C[2] = real(dot(ϕ_l, ϕ_1))
    #
    for m in 3:M
        ϕ_2 = 2 * H_scaled * ϕ_1 - ϕ_0
        C[m] = real(dot(ϕ_l, ϕ_2))
        ϕ_0[:] .= ϕ_1[:]
        ϕ_1[:] .= ϕ_2[:]
    end
    return C
end

function sum_chebyshev(M, C, E_scaled, E_max)
    g = jackson_damping(M) # get the damping factor
    C .*= g                 # apply the damping factor
    Ne = length(E_scaled)     # number of energy points
    T0 = ones(Ne)
    T1 = copy(E_scaled)
    T2 = zeros(Ne)
    s = zeros(Ne)
    for i in 1:Ne
        s[i] = C[2] * T1[i]
    end
    for m in range(3, M)
        for i in 1:Ne
            T2[i] = 2.0 * E_scaled[i] * T1[i] - T0[i]
        end
        T0[:] .= T1[:]
        T1[:] .= T2[:]
        for i in 1:Ne
            s[i] += C[m] * T2[i]
        end
    end
    s .*= 2.0
    s .+= C[1]
    for i in 1:Ne
        s[i] *= 2.0 / (pi * sqrt(1.0 - E_scaled[i]^2) * E_max)
    end
    return s
end

function calc_dos(M, E_max, E_scaled, H_scaled, phi)
    C = calc_moments(M, H_scaled, phi, phi)
    DOS = sum_chebyshev(M, C, E_scaled, E_max)
    return DOS
end

function evolve(H_scaled, dt_scaled, sign, phi_i)
    phi_0 = phi_i
    phi_1 = H_scaled * phi_i
    j0 = besselj(0, dt_scaled)
    j1 = besselj(1, dt_scaled)
    phi_o = j0 * phi_0 + 2.0 * (-1im * sign) * j1 * phi_1
    m = 2
    while true
        jm = besselj(m, dt_scaled)
        if abs(jm) < 1.0e-15
            break
        end
        phi_2 = 2.0 * H_scaled * phi_1 - phi_0
        phi_o += 2.0 * (-im*sign)^m * jm * phi_2
        phi_0 = phi_1
        phi_1 = phi_2
        m += 1
    end
    return phi_o
end

function calc_vac(M, E_max, dt_scaled, E_scaled, H_scaled, V, ϕ, DOS)
    Ne = length(E_scaled)  # number of energy points
    Nt = length(dt_scaled) # number of time steps
    dt = dt_scaled / E_max
    phi_left = ϕ
    phi_right = V * ϕ
    DOS_times_VAC_old = zeros(Ne)
    DOS_times_VAC_new = zeros(Ne)
    VAC = zeros( Nt, Ne )
    sigma_from_VAC = zeros((Nt, Ne))
    for nt in 1:Nt
        C = calc_moments(M, H_scaled, V * phi_left, phi_right)
        DOS_times_VAC_new = sum_chebyshev(M, C, E_scaled, E_max)
        VAC[nt, :] = DOS_times_VAC_new ./ DOS
        if nt > 1
            tmp = dt[nt-1] * (DOS_times_VAC_old + DOS_times_VAC_new) * 0.5
            sigma_from_VAC[nt, :] = sigma_from_VAC[nt - 1, :] + tmp
        end
        DOS_times_VAC_old[:] = DOS_times_VAC_new
        phi_left[:] = evolve(H_scaled, dt_scaled[nt], -1, phi_left)
        phi_right[:] = evolve(H_scaled, dt_scaled[nt], -1, phi_right)
    end
    sigma_from_VAC *= 2.0 * pi # from e^2/hbar to e^2/h
    return VAC, sigma_from_VAC
end


function evolvex(H_scaled, V_scaled, dt_scaled, phi_i)
    phi_0 = copy(phi_i)
    phix_0 = zeros(size(phi_i)) #phi_i * 0.0
    phi_1 = H_scaled * phi_0
    phix_1 = im * V_scaled * phi_0
    phi_o = -2im * besselj(1, dt_scaled) * phix_1
    m = 2
    while true
        jm = besselj(m, dt_scaled)
        if abs(jm) < 1.0e-15
            break
        end
        phi_2 = 2.0 * H_scaled * phi_1 - phi_0
        phix_2 = 2im * V_scaled * phi_1
        phix_2 += 2.0 * H_scaled * phix_1 - phix_0
        phi_o += 2.0 * (-1im)^m * jm * phix_2
        phi_0 = phi_1 
        phi_1 = phi_2
        phix_0 = phix_1
        phix_1 = phix_2
        m += 1
    end
    return phi_o
end


function calc_msd(M, E_max, dt_scaled, E_scaled, H_scaled, V_scaled, phi, DOS)
    Ne = length(E_scaled)  # number of energy points
    Nt = length(dt_scaled) # number of time steps
    dt = dt_scaled / E_max
    MSD = zeros((Nt, Ne))
    sigma_from_MSD = zeros((Nt, Ne))
    phix = zeros(size(phi))
    DOS_times_MSD_old = zeros(Ne)
    DOS_times_MSD_new = zeros(Ne)
    for nt in 1:Nt
        phix = evolve(H_scaled, dt_scaled[nt], 1, phix)
        phix += evolvex(H_scaled, V_scaled, dt_scaled[nt], phi)
        phi = evolve(H_scaled, dt_scaled[nt], 1, phi)
        C = calc_moments(M, H_scaled, phix, phix)
        DOS_times_MSD_new = sum_chebyshev(M, C, E_scaled, E_max);
        MSD[nt, :] = DOS_times_MSD_new ./ DOS
        sigma_from_MSD[nt, :] = (DOS_times_MSD_new - DOS_times_MSD_old) / dt[nt] * 0.5
        DOS_times_MSD_old = DOS_times_MSD_new
    end
    sigma_from_MSD *= 2.0 * pi # from e^2/hbar to e^2/h
    return MSD, sigma_from_MSD
end

