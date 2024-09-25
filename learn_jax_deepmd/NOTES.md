Training and validation separation is done separately from training code.

Data format can be read directly using Numpy.

Load data:
```python
train_data = data.DPDataset(train_paths, labels, {"atomic_sel": atomic_sel})
```

What is `atomic_sel`? This is only used for per-atom property prediction?


