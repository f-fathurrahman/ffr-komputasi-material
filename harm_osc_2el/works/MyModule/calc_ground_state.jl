function calc_ground_state!(state; max_iter=1000, tol=1e-10, verbose = 0)
    
    @info "--------------------------------"
    @info "Calculating ground state"
    @info "typeof(state) = $(typeof(state))"
    @info "--------------------------------"
    
    old_E = energy(state)
    delta_E = 0
    for i in 1:max_iter
        update!(state)
        #
        new_E = energy(state)
        delta_E = old_E - new_E
        #
        if abs(delta_E) < tol
            if verbose != 0
                @printf("Converged !!!")
            end
            return
        end
        if verbose == 1
            @printf("Iteration %5d   E = %18.10f    ΔE = %15.5e\n", i, new_E, delta_E)
        end
        old_E = new_E
    end
    println("Did not converge after $(max_iter) iterations. Final energy change was $(delta_E)")
    return
end