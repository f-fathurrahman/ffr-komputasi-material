from prepare_args_01 import args
import my_jrystal as jr


def main_debug():
    print(args)
    config = jr.config.get_config(args.config)
    # Calculate energy
    jr.calc.energy_normcons(config)
    #jr.calc.energy_all_electrons(config)


# XXX: We must write it like this due to some multiprocessing stuffs
if __name__ == "__main__":
    main_debug()
