import warnings
warnings.filterwarnings("ignore")
from my_mace.cli.run_train import main as mace_run_train_main
import sys
import logging

def train_mace(config_file_path):
    logging.getLogger().handlers.clear()
    # modify program args
    sys.argv = ["program", "--config", config_file_path]
    mace_run_train_main()

# Do the training
train_mace("config_v01.yml")
