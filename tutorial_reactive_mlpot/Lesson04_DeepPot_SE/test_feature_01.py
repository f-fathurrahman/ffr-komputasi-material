import math
from typing import Sequence, Tuple

import numpy as np

import torch
from torch import Tensor
from torch import nn

from load_data import load_data
from models import *

from local_environment import local_environment

from  feature import Feature


coord, atom_types, grad = load_data()

descriptor = Feature(4, neuron=[25, 50], axis_neuron=4)

# only one instance of atom_types are used
res = descriptor(coord[1:4], atom_types[0])
print(res)

# atom_types does not change

print("Pass here ...")
