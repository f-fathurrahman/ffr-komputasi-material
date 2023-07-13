import Random

data_keys = (energy_key = "energy", force_key = "forces")

function gen_dat()
   sw = JuLIP.Potentials.StillingerWeber() 
   at = JuLIP.Utils.rattle!(
      JuLIP.Build.bulk(:Si, cubic=true) * rand(2:3), 0.3
   )
   JuLIP.set_data!(at, data_keys.energy_key, JuLIP.energy(sw,at))
   JuLIP.set_data!(at, data_keys.force_key, JuLIP.forces(sw,at))
   return at
end

Random.seed!(0)
train = [gen_dat() for _=1:20];