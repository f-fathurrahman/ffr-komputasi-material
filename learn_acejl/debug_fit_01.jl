using MyACEpotentials

raw_data1 = read_extxyz("datasets/TiAl_tutorial.xyz");
model1 = acemodel(
    elements = [:Ti, :Al],
    order = 3, totaldegree = 12,
    Eref = [:Ti => -1586.0195, :Al => -105.5954]
)

#=
raw_data2 = read_extxyz("../pyxtal_ff/TEMP_DATASET/Ni_fcc_4atoms_300K_v1.xyz");
model2 = acemodel(
    elements = [:Ni],
    order = 3, totaldegree = 12,
    Eref = [:Ni => -1000.0]
);
=#

# from default_weights ?
weights = Dict(
    "default" => Dict("E"=>1.0, "F"=>1.0, "V"=>1.0)
);

data1 = map( raw_data1 ) do data_point
    MyACEpotentials._apply_weight(
        data_point;
        energy_key = "energy", 
        force_key = "force", 
        virial_key = "virial", 
        weights = weights, 
        v_ref = model1.Vref
    )
end

#=
data2 = map( raw_data2 ) do data_point
    MyACEpotentials._apply_weight(
        data_point;
        energy_key = "energy", 
        force_key = "forces", 
        virial_key = "stress", 
        weights = weights, 
        v_ref = model2.Vref
    )
end
=#


