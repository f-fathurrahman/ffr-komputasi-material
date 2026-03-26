# DMC: Diffusion Monte Carlo Programs

This directory contains various DMC programs and related utilities.

# Main Programs
## Python programs are short test programs used in the book,
	- `expH_sympy.py`
    - `harm_osc_density_matrix_tests.py`
    - `particle_in_a_box_no_importance_sampling.py`
	- `sympy_drift.py`
	- `two_bosons_in_a_box_no_importance_sampling.py`
    - `H2_bonding_antibonding.py`
	- `particle_in_a_box_drift_trajectories.py`
	- `particle_in_a_box.py`
    - `triple_sympy.py`

## Julia programs 
  - `two_fermions_in_a_box_no_importance_sampling.jl`
      - Standalone program used in the book
  - `dmc.heatom.jl` 
	  - Standalone program basic DMC for He atom 
  - `Atom_Slater_Jastrow.jl` 
	- DMC using Slater-Jastrow trial wave function
		- Be atom using STO basis set, run \		
		 `julia Atom_Slater_Jastrow.jl atom=Be basis_set=STO`		
		- Input data is in file `atom_data/Be` and the orbital parameters in file `atom_data/Be_STO`
			- Only H, He, Be, B and Ne 
			- orbital parameter file contains STO coefficients and exponents
			for each orbital (1s, 2s, 2p ...),
			Jastrow factor parameters, and finally 
			the linear combination coefficients (`lc coefficients`)\
			where `coeff` is the linear coefficient of the orbitals listed on the next line,
			for example\
			`orbitals = 10 10 21 21`\
			are the `lm` quantum numbers of 1s, 1s, 2p, and 2p orbitals.\
			If the input file is a result of optimization, it contains also `raw output of wf_params` (optimal wave function 
			parameters) and the estimated energy computed in optimization.

        - Similarly for double-zeta STO basis (basis_set=STO_Dzeta) 
		and triple-zeta basis (basis_set=STO_Tzeta) 
	
	- `Atom_Slater_Jastrow_optimization.jl` 
		- optimization of basis set parameters using initial parameters in `atom_data/`
			- To optimize B atom using triple-zeta STO basis, run\
				`julia Atom_Slater_Jastrow_optimization.jl atom=B basis_set=Tzeta_STO`
				- raw output to file Eopt_B_STO_Tzeta
				- optimized input is stored to file `atom_data/B_Tzeta_STO_opt`, which can be copied to 
				`atom_data/B_Tzeta_STO` and used as input.
		 
	- `H2_Kolos_Roothaan.jl`
		- DMC on H<sub>2</sub> molecule using the Kolos-Roothaan trial wave function.
			- Optimize parameter for fixed proton-proton distance Req=1.2 using\
			`julia H2_Kolos_Roothaan.jl optimize=1 Req=1.2`\
			copy the parameters to the code (which already has the data used in the book)
		- Run DMC using $\tau$=0.001,\
			`julia H2_Kolos_Roothaan.jl  tau=0.001 Req=1.2`

# Supporting Julia Modules
- `Utilities.jl`
- `QMC_Statistics.jl` function to collect data samples, block them, and give results. 
- `HydrogenicOrbitals.jl`, `STOOrbitals.jl`, `SlaterDeterminants.jl` orbitals and building Slater determinants
- `LinearOptimization.jl` programs for linear optimization of trial wave function parameters
- `VMCstep.jl` a few variants of performing a single VMC step
- `Common.jl` a collection of common utility programs used in Slater-Jastrow codes
