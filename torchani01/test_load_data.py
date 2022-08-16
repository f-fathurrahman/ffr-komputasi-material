import os
import sys
import torch
import my_torchani

torch.manual_seed(1234) # for reproducibility

# ani_gdb_s01.h5  ani_gdb_s02.h5  ani_gdb_s03.h5  ani_gdb_s04.h5
#dspath = "/home/efefer/WORKS/TORCHANI/torchani-master/dataset/ani1-up_to_gdb4/ani_gdb_s01.h5"
dspath = "/home/efefer/WORKS/TORCHANI/torchani-master/dataset/ani-1x/sample.h5"
dataset = my_torchani.data.load(dspath, additional_properties=('forces',))

print("type(dataset) = ", type(dataset)) # TransformableIterable

DATASET_LIST = []
for i,v in enumerate(dataset):
    DATASET_LIST.append(v)
print("Length DATASET_LIST = ", len(DATASET_LIST))
print("Sizeof DATASET_LIST = %f MB" % (sys.getsizeof(DATASET_LIST)/1024/1024))

# Now DATASET_LIST contains list of dict
print("Keys of dict in DATASET_LIST:")
print(DATASET_LIST[0].keys())
