function compute_ground_state!(state; max_iter = 1000, tol = 1e-10, verbose = 0)
    old_E = energy(state)
    delta_E = 0
    for i in 1:max_iter
        update!(state)
        
        new_E = energy(state)
        delta_E = old_E - new_E

        if abs(delta_E) < tol
            if verbose != 0
                println("Iteration $i: E = $new_E")
            end

            return state
        end
        if verbose == 1 || (verbose == 2 && (i <= 3 || i >= max_iter - 3))
            println("Iteration $i: E = $new_E")
        end
        if verbose == 3
            print("\rIteration $i: E = $new_E                           ")
        end
        old_E = new_E
    end
    println("Did not converge after $(max_iter) iterations. Final energy change was $(delta_E)")
    return state
end