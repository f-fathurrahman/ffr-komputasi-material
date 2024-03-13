using DFTK
using Plots
using Unitful
using UnitfulAtomic

a = 3.01794  # bohr
b = 5.22722  # bohr
c = 9.77362  # bohr
lattice = [[-a -a  0]; [-b  b  0]; [0   0 -c]]
Mg = ElementPsp(:Mg; psp=load_psp("hgh/pbe/Mg-q2"))
atoms     = [Mg, Mg]
positions = [[2/3, 1/3, 1/4], [1/3, 2/3, 3/4]];

kspacing = 0.945 / u"angstrom"        # Minimal spacing of k-points,
#                                      in units of wavevectors (inverse Bohrs)
Ecut = 5                              # Kinetic energy cutoff in Hartree
temperature = 0.01                    # Smearing temperature in Hartree
smearing = DFTK.Smearing.FermiDirac() # Smearing method
#                                      also supported: Gaussian,
#                                      MarzariVanderbilt,
#                                      and MethfesselPaxton(order)

model = model_DFT(lattice, atoms, positions, [:gga_x_pbe, :gga_c_pbe];
                  temperature, smearing)
kgrid = kgrid_from_maximal_spacing(lattice, kspacing)
basis = PlaneWaveBasis(model; Ecut, kgrid);

scfres = self_consistent_field(basis, damping=0.8, mixing=KerkerMixing());


