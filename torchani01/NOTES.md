
The data is stored in h5 a file. This file can be loaded my calling:

```python
dspath = "/home/efefer/WORKS/TORCHANI/torchani-master/dataset/ani1-up_to_gdb4/ani_gdb_s01.h5"
dataset = torchani.data.load(dspath)
```


Batching:
```python
batch_size = 2560
training = training.collate(batch_size).cache()
validation = validation.collate(batch_size).cache()
```

Padding of coordinate and species.
For coordinate zero padding is used. For species invalid species index is used, in this case
it is -1.
```py
In [35]: DATASET_TRAINING[0]["coordinates"][2]                                                                                         
Out[35]: 
tensor([[-9.5886e-19, -2.8211e-03,  1.0855e-01],
        [ 1.0322e-17,  8.4906e-01, -4.1226e-01],
        [ 4.8957e-18, -8.0428e-01, -3.7733e-01],
        [ 0.0000e+00,  0.0000e+00,  0.0000e+00],
        [ 0.0000e+00,  0.0000e+00,  0.0000e+00]])

In [36]: DATASET_TRAINING[0]["species"][2]                                                                                             
Out[36]: tensor([ 3,  0,  0, -1, -1])
```