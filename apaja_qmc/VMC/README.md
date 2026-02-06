# VMC: Variational Monte Carlo Programs

This directory contains mostly VMC programs.


# Main Programs

- `vmc.heatom_direct_vs_correlated_sampling.jl` 
  - Demonstrates the difference of direct and correlated sampling 
  - The package "Plots" loads very slowly, so if you run it several times, use the Julia REPL

- `vmc.heatom_optimization.jl` optimizes a He atom trial wave function

- `vmc.heatom.jl` computes the energy of a He atom using the trial wave function defined in Model_Heatom.jl

# Helpers

- `Common.jl` defines the structure VMC_Params for acceptance and step size in VMC (and many other things for DMC)

- `Utilities.jl` defines output format of QMC results (among other things)

- `QMC_Statistics.jl` functions to collect data and evaluate error estimates
