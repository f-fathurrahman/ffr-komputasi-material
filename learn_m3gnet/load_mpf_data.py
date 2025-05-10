import random
import numpy as np

import pickle as pk
import pandas as pd
import pymatgen

# From: https://gitlab.brindisi.enea.it/claudio.ronchetti/ai4mat/-/blob/c8864a90a12b8fb32f95285e27db531d12c3f856/ml/alignn/data.py
def get_id_train_val_test(
    total_size=1000,
    split_seed=123,
    train_ratio=None,
    val_ratio=0.1,
    test_ratio=0.1,
    n_train=None,
    n_test=None,
    n_val=None,
    keep_data_order=False,
):
    """Get train, val, test IDs."""
    if (
        train_ratio is None
        and val_ratio is not None
        and test_ratio is not None
    ):
        if train_ratio is None:
            assert val_ratio + test_ratio < 1
            train_ratio = 1 - val_ratio - test_ratio
            print("Using rest of the dataset except the test and val sets.")
        else:
            assert train_ratio + val_ratio + test_ratio <= 1
    # indices = list(range(total_size))
    if n_train is None:
        n_train = int(train_ratio * total_size)
    if n_test is None:
        n_test = int(test_ratio * total_size)
    if n_val is None:
        n_val = int(val_ratio * total_size)
    ids = list(np.arange(total_size))
    if not keep_data_order:
        random.seed(split_seed)
        random.shuffle(ids)
    # np.random.shuffle(ids)
    if n_train + n_val + n_test > total_size:
        raise ValueError(
            "Check total number of samples.",
            n_train + n_val + n_test,
            ">",
            total_size,
        )

    # shuffle consistently with https://github.com/txie-93/cgcnn/data.py
    # i.e. shuffle the index in place with standard library random.shuffle
    # first obtain only valid indices

    # test_size = round(N * 0.2)

    # full train/val test split
    # ids = ids[::-1]
    id_train = ids[:n_train]
    id_val = ids[-(n_val + n_test) : -n_test]  # noqa:E203
    id_test = ids[-n_test:]
    return id_train, id_val, id_test



print('loading the MPF dataset 2021')
with open('MPF_data/block_0.p', 'rb') as f:
    data = pk.load(f)

#with open('MPF_data/block_1.p', 'rb') as f:
#    data2 = pk.load(f)

print('MPF dataset 2021 loaded')

"""
#data.update(data2)
df = pd.DataFrame.from_dict(data)

id_train, id_val, id_test = get_id_train_val_test(
    total_size=len(data),
    split_seed=42,
    train_ratio=0.90,
    val_ratio=0.05,
    test_ratio=0.05,
    keep_data_order=False,
)

dataset_train = {}
dataset_val = {}
dataset_test = {}
cnt = 0
for idx, item in df.items():
    # import pdb; pdb.set_trace()
    if cnt in id_train:
        for iid in range(len(item['energy'])):
            dataset_train.append({
                "atoms" : item['structure'][iid],
                "energy" : item['energy'][iid] / len(item['force'][iid]),
                "force": np.array(item['force'][iid])
            })
    elif cnt in id_val:
        for iid in range(len(item['energy'])):
            dataset_val.append({
                "atoms" : item['structure'][iid],
                "energy" : item['energy'][iid] / len(item['force'][iid]),
                "force": np.array(item['force'][iid])})
    elif cnt in id_test:
        for iid in range(len(item['energy'])):
            dataset_test.append({
                "atoms": item['structure'][iid],
                "energy": item['energy'][iid] / len(item['force'][iid]),
                "force": np.array(item['force'][iid])
            })
    cnt += 1

print('using %d samples to train, %d samples to evaluate, and %d samples to test ' % (len(dataset_train), len(dataset_val), len(dataset_test)))
"""


