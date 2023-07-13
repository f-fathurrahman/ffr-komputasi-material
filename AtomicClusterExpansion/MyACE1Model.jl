mutable struct MyACE1Model 
   basis 
   params 
   Vref  
   potential
   meta::Dict
end


function my_acemodel(;  kwargs...)
   Eref = get(kwargs, :Eref, nothing)
   Vref = get(kwargs, :Vref, nothing)

   println("Eref = ", Eref)
   println("Eref = ", Eref)

   # construct the basis 
   basis = ACE1x.ace_basis(; kwargs...)

   # Not executed
   if Vref == nothing && Eref != nothing 
      Vref = JuLIP.OneBody(Eref...)
   end

   # construct a model without parameters and without an evaluator 
   model = MyACE1Model(basis, nothing, Vref, nothing, Dict())

   # set some random parameters -> this will also generate the evaluator 
   params = randn(length(basis))
   params = params ./ (1:length(basis)).^4
   _set_params!(model, params)
   return model 
end


_sumip(pot1, pot2) = 
      JuLIP.MLIPs.SumIP([pot1, pot2])

_sumip(pot1::JuLIP.MLIPs.SumIP, pot2) = 
      JuLIP.MLIPs.SumIP([pot1.components..., pot2])
                                      
function _set_params!(model, params)
   model.params = params
   model.potential = JuLIP.MLIPs.combine(model.basis, model.params)
   if model.Vref != nothing 
      model.potential = _sumip(model.potential, model.Vref)
   end
   return model # XXX why need this? 
end



function my_acefit!(
   model::MyACE1Model, raw_data;
   solver = ACEfit.BLR(),
   weights = ACE1pack.default_weights(),
   energy_key = "energy", 
   force_key = "force", 
   virial_key = "virial", 
   smoothness = 2, 
   prior = nothing, 
   repulsion_restraint = false, 
   restraint_weight = 0.01, 
   export_lammps = nothing, 
   export_json = nothing, 
   verbose=true
)

   data = [ ACE1pack.AtomsData(at; energy_key = energy_key, force_key=force_key, 
                          virial_key = virial_key, weights = weights, 
                          v_ref = model.Vref) 
            for at in raw_data ] 

   if verbose
      ACE1pack.assess_dataset(data)
   end 

   if repulsion_restraint 
      append!(data, _rep_dimer_data(model, weight = restraint_weight))
   end
                  
   P = ACE1pack._make_prior(model.basis, smoothness, prior)
   println("Pass here")

   A, Y, W = ACEfit.assemble(data, model.basis)
   Ap = Diagonal(W) * (A / P) 
   Y = W .* Y
   
   result = ACEfit.solve(solver, Ap, Y)
   coeffs = P \ result["C"]
   ACE1x._set_params!(model, coeffs)

   return model 
end



function my_linear_errors(
    raw_data::AbstractVector{<: JuLIP.Atoms}, model::MyACE1Model; 
    energy_key = "energy", 
    force_key = "force", 
    virial_key = "virial", 
    weights = ACE1pack.default_weights()
)
    #
    Vref = model.Vref                       
    data = [ ACE1pack.AtomsData(at; energy_key = energy_key, force_key=force_key, 
                          virial_key = virial_key, weights = weights, 
                          v_ref = model.Vref) 
            for at in raw_data ] 
   return my_linear_errors(data, model.potential)
end


function my_linear_errors(data, model; group_key="config_type")

   mae = Dict("E"=>0.0, "F"=>0.0, "V"=>0.0)
   rmse = Dict("E"=>0.0, "F"=>0.0, "V"=>0.0)
   num = Dict("E"=>0, "F"=>0, "V"=>0)

   config_types = []
   config_mae = Dict{String,Any}()
   config_rmse = Dict{String,Any}()
   config_num = Dict{String,Any}()

   for d in data

       c_t = ACE1pack.group_type(d; group_key)
       if !(c_t in config_types)
          push!(config_types, c_t)
          merge!(config_rmse, Dict(c_t=>Dict("E"=>0.0, "F"=>0.0, "V"=>0.0)))
          merge!(config_mae, Dict(c_t=>Dict("E"=>0.0, "F"=>0.0, "V"=>0.0)))
          merge!(config_num, Dict(c_t=>Dict("E"=>0, "F"=>0, "V"=>0)))
       end

       # energy errors
       if !isnothing(d.energy_key)
           estim = JuLIP.energy(model, d.atoms) / length(d.atoms)
           exact = d.atoms.data[d.energy_key].data / length(d.atoms)
           mae["E"] += abs(estim-exact)
           rmse["E"] += (estim-exact)^2
           num["E"] += 1
           config_mae[c_t]["E"] += abs(estim-exact)
           config_rmse[c_t]["E"] += (estim-exact)^2
           config_num[c_t]["E"] += 1
       end

       # force errors
       if !isnothing(d.force_key)
           estim = JuLIP.mat(JuLIP.forces(model, d.atoms))
           exact = JuLIP.mat(d.atoms.data[d.force_key].data)
           mae["F"] += sum(abs.(estim-exact))
           rmse["F"] += sum((estim-exact).^2)
           num["F"] += 3*length(d.atoms)
           config_mae[c_t]["F"] += sum(abs.(estim-exact))
           config_rmse[c_t]["F"] += sum((estim-exact).^2)
           config_num[c_t]["F"] += 3*length(d.atoms)
       end

       # virial errors
       if !isnothing(d.virial_key)
           estim = virial(model, d.atoms)[SVector(1,5,9,6,3,2)] ./ length(d.atoms)
           exact = d.atoms.data[d.virial_key].data[SVector(1,5,9,6,3,2)] ./ length(d.atoms)
           mae["V"] += sum(abs.(estim-exact))
           rmse["V"] += sum((estim-exact).^2)
           num["V"] += 6
           config_mae[c_t]["V"] += sum(abs.(estim-exact))
           config_rmse[c_t]["V"] += sum((estim-exact).^2)
           config_num[c_t]["V"] += 6
       end
    end

    # finalize errors
    for (k,n) in num
        (n==0) && continue
        rmse[k] = sqrt(rmse[k] / n)
        mae[k] = mae[k] / n
    end
    errors = Dict("mae"=>mae, "rmse"=>rmse)

    # finalize config errors
    for c_t in config_types
        for (k,c_n) in config_num[c_t]
            (c_n==0) && continue
            config_rmse[c_t][k] = sqrt(config_rmse[c_t][k] / c_n)
            config_mae[c_t][k] = config_mae[c_t][k] / c_n
        end
    end
    config_errors = Dict("mae"=>config_mae, "rmse"=>config_rmse)

    # merge errors into config_errors and return
    push!(config_types, "set")
    merge!(config_errors["mae"], Dict("set"=>mae))
    merge!(config_errors["rmse"], Dict("set"=>rmse))

    @info "RMSE Table"
    header = ["Type", "E [meV]", "F [eV/A]", "V [meV]"]
    table = hcat(
        config_types,
        [1000*config_errors["rmse"][c_t]["E"] for c_t in config_types],
        [config_errors["rmse"][c_t]["F"] for c_t in config_types],
        [1000*config_errors["rmse"][c_t]["V"] for c_t in config_types],
    )
    PrettyTables.pretty_table(
        table; header=header,
        body_hlines=[length(config_types)-1],
        formatters=PrettyTables.ft_printf("%5.3f"),
        crop = :horizontal)

    @info "MAE Table"
    header = ["Type", "E [meV]", "F [eV/A]", "V [meV]"]
    table = hcat(
        config_types,
        [1000*config_errors["mae"][c_t]["E"] for c_t in config_types],
        [config_errors["mae"][c_t]["F"] for c_t in config_types],
        [1000*config_errors["mae"][c_t]["V"] for c_t in config_types],
    )
    PrettyTables.pretty_table(
        table; header=header,
        body_hlines=[length(config_types)-1],
        formatters=PrettyTables.ft_printf("%5.3f"),
        crop = :horizontal)

    return config_errors
end