using MyACEpotentials

#raw_data = read_extxyz("../pyxtal_ff/TEMP_DATASET/Ni_fcc_4atoms_300K_v1.xyz");
raw_data = read_extxyz("../pyxtal_ff/WORKS_Ni_01/DATASET_26062025/DATA_Ni_26_06_2025.xyz");
model = acemodel(
    elements = [:Ni],
    order = 3, totaldegree = 12,
    #Eref = [:Ni => -1000.0]
);

acefit!(model, raw_data,
    energy_key = "energy", 
    force_key = "forces", 
    virial_key = "stress", 
);
errors = MyACEpotentials.linear_errors(raw_data, model)
