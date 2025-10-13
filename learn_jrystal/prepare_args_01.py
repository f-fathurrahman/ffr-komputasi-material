import argparse

parser = argparse.ArgumentParser(
  prog='Jrystal', description='Command for Jrystal package.'
)

parser.add_argument(
  "-m",
  "--mode",
  choices=["energy", "band"],
  default='energy',
  help="Set the computation mode. For total enrgy minimization, please use "
  "\'energy\'. For band structure calculation, please use \'band\'. "
)

# Modify this
parser.add_argument(
  "-c",
  "--config",
  default='Al_config.yaml',
  help="Set the configuration file path."
)

parser.add_argument(
  "-l",
  "--load",
  help=(
    "Load pickled output from energy calculation for band structure "
    "calculation."
  )
)

args = parser.parse_args()
