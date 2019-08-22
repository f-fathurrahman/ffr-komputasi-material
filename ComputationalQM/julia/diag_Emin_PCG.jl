#=
function diag_Emin_PCG( Ham::Hamiltonian, X0::Array{ComplexF64,2};
                        tol=1e-5, NiterMax=100, verbose=false,
                        verbose_last=false, Nstates_conv=0,
                        tol_ebands=1e-4,
                        α_t=3e-5, I_CG_BETA = 2 )

    X = copy(X0)
    evals = diag_Emin_PCG!( Ham, X, tol=tol, NiterMax=NiterMax, verbose=verbose,
                            verbose_last=verbose_last, Nstates_conv=Nstates_conv,
                            tol_ebands=tol_ebands,
                            α_t=α_t, I_CG_BETA = I_CG_BETA )                          
    return evals, X
end


function diag_Emin_PCG!( Ham::Hamiltonian, psiks::BlochWavefunc;
                         tol=1e-5, NiterMax=100, verbose=false,
                         verbose_last=false, Nstates_conv=0,
                         tol_ebands=1e-4,
                         α_t=3e-5, I_CG_BETA = 2 )
    
    pw = Ham.pw
    electrons = Ham.electrons
    Nkpt = pw.gvecw.kpoints.Nkpt
    Nspin = electrons.Nspin
    Nkspin = Nspin*Nkpt
    Ngw = pw.gvecw.Ngw
    Nstates = electrons.Nstates

    evals = zeros(Float64,Nstates,Nkspin)

    for ispin = 1:Nspin
    for ik = 1:Nkpt
        Ham.ik = ik
        Ham.ispin = ispin
        ikspin = ik + (ispin - 1)*Nkpt
        #
        evals[:,ikspin] =
        diag_Emin_PCG!( Ham, psiks[ikspin], tol=tol, NiterMax=NiterMax, verbose=verbose,
                       verbose_last=verbose_last, Nstates_conv=Nstates_conv,
                       tol_ebands=tol_ebands,
                       α_t=α_t, I_CG_BETA = I_CG_BETA )
        #
    end
    end

    return evals

end
=#


# used by diag_Emin_PCG
function calc_grad_evals( Ham, psi::Array{Float64,2} )
    Nbasis  = size(psi,1)
    Nstates = size(psi,2)
    #
    grad = zeros(Float64,Nbasis,Nstates)
    #
    H_psi = Ham*psi
    for ist in 1:Nstates
        grad[:,ist] = H_psi[:,ist]
        for jst in 1:Nstates
            grad[:,ist] = grad[:,ist] - dot( psi[:,jst], H_psi[:,ist] ) * psi[:,jst]
        end
    end
    return grad
end


function diag_Emin_PCG!( Ham, X::Array{Float64,2};
                         tol=1e-5,
                         NiterMax=100,
                         verbose=false,
                         verbose_last=false,
                         Nstates_conv=0,
                         tol_ebands=1e-4,
                         α_t=3e-5,
                         i_cg_beta=2 )

    Nbasis  = size(X,1)
    Nstates = size(X,2)

    if Nstates_conv == 0
        Nstates_conv = Nstates
    end

    #
    # Variabls for PCG
    #
    d = zeros(Float64,Nbasis,Nstates)
    g_old = copy(d)
    d_old = copy(d)
    Kg = copy(d)
    Kg_old = copy(d)
    Xc = copy(d)
    gt = copy(d)
    
    β = 0.0
    
    Ebands_old = 0.0

    Hr = Symmetric( X' * (Ham*X) )

    evals = eigvals(Hr)
    Ebands = sum(evals)
    
    evals_old = copy(evals)
    devals = ones(Nstates)

    IS_CONVERGED = false

    for iter = 1:NiterMax

        g = calc_grad_evals( Ham, X)
        #Kg = Kprec( Ham.ik, pw, g )
        Kg = copy(g)

        if iter != 1
            if i_cg_beta == 1
                β = real(sum(conj(g).*Kg))/real(sum(conj(g_old).*Kg_old))
            elseif i_cg_beta == 2
                β = real(sum(conj(g-g_old).*Kg))/real(sum(conj(g_old).*Kg_old))
            elseif i_cg_beta == 3
                β = real(sum(conj(g-g_old).*Kg))/real(sum(conj(g-g_old).*d))
            else
                β = real(sum(conj(g).*Kg))/real(sum((g-g_old).*conj(d_old)))
            end
        end
        if β < 0.0
            β = 0.0
        end

        d = -Kg + β * d_old

        Xc = ortho_gram_schmidt( X + α_t*d )
        gt = calc_grad_evals( Ham, Xc )

        denum = real(sum(conj(g-gt).*d))
        if denum != 0.0
            α = abs( α_t*real(sum(conj(g).*d))/denum )
        else
            α = 0.0
        end

        # Update wavefunction
        X = X + α*d
        ortho_gram_schmidt!(X)

        Hr = Symmetric( X' * ( Ham*X ) )
        evals = eigvals(Hr)
        Ebands = sum(evals)

        devals = abs.( evals - evals_old )
        evals_old = copy(evals)
        
        nconv = length( findall( devals .< tol ) )

        diffE = abs(Ebands-Ebands_old)

        if verbose
            @printf("CG step %8d = %18.10f %10.7e\n", iter, Ebands, diffE)
            for ist = 1:Nstates
                @printf("evals[%3d] = %18.10f, devals = %18.10e\n", ist, evals[ist], devals[ist] )
            end
            @printf("iter %d nconv = %d\n", iter, nconv)
        end
        if nconv >= Nstates_conv
            IS_CONVERGED = true
            if verbose || verbose_last
                @printf("Convergence is achieved based on nconv\n")
            end
            break
        end
        if diffE <= tol_ebands*Nstates
            IS_CONVERGED = true
            if verbose || verbose_last
                @printf("Convergence is achieved based on tol_ebands*Nstates\n")
            end
            break
        end

        g_old = copy(g)
        d_old = copy(d)
        Kg_old = copy(Kg)
        Ebands_old = Ebands
    end

    if !IS_CONVERGED
        @printf("\nWARNING: diag_Emin_PCG is not converged after %d iterations\n", NiterMax)
    end

    ortho_gram_schmidt!(X)
    Hr = Symmetric( X' * (Ham*X) )
    evals, evecs = eigen(Hr)
    # Sort
    idx_sorted = sortperm(evals)
    evals = evals[idx_sorted]
    X[:,:] = X*evecs[:,idx_sorted]

    if verbose_last || verbose
        @printf("\nEigenvalues from diag_Emin_PCG:\n\n")
        for ist = 1:Nstates
            @printf("evals[%3d] = %18.10f devals = %18.10e\n", ist, evals[ist], devals[ist] )
        end
    end

    return evals
end
