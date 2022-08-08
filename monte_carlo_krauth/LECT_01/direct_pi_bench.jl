using Printf
using BenchmarkTools
using Distributions

function direct_pi_v1(N)
    Nhits = 0
    for i in 1:N
        x = 2*( rand() - 0.5 )
        y = 2*( rand() - 0.5 )
        if x^2 + y^2 < 1.0
            Nhits += 1
        end
    end
    return
end

function direct_pi_v2(N)
    Nhits = 0
    rnd = Distributions.Uniform(-1.0, 1.0)
    for i in 1:N
        x = rand(rnd)
        y = rand(rnd)
        if x^2 + y^2 < 1.0
            Nhits += 1
        end
    end
    return
end

@btime direct_pi_v1(5000)

@btime direct_pi_v2(5000)

