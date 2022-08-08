using Printf

function direct_pi(N)
    Nhits = 0
    for i in 1:N
        x = 2*( rand() - 0.5 )
        y = 2*( rand() - 0.5 )
        if x^2 + y^2 < 1.0
            Nhits += 1
        end
    end
    return Nhits
end


function main()
    Nruns = 1000
    Ntrials = 1_000_000
    for irun in 1:Nruns
        Nhits = direct_pi(Ntrials)
        pi_approx = 4.0 * Nhits / Ntrials
        err = abs(pi - pi_approx)
        @printf("%5d  %18.10f  %10.5e\n", irun, pi_approx, err)
    end
end

main()
