from prepare_args_01 import args
print(args)

import my_jrystal as jr
config = jr.config.get_config(args.config)

# Calculate energy
jr.calc.energy_normcons(config)
#jr.calc.energy_all_electrons(config)