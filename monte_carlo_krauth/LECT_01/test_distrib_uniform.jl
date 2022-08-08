using Distributions

rnd = Distributions.Uniform(-1.0, 1.0)

Nsamples = 1_000_000
x = rand(rnd, Nsamples)
y = rand(rnd, Nsamples)
in_circle = (x.^2 + y.^2) .< 1.0
Nhits = sum(in_circle)

pi_approx = 4*Nhits/Nsamples
println("pi_approx = ", pi_approx)
