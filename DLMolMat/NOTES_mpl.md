```python
import matplotlib.pyplot as plt
import matplotlib
matplotlib.style.use("seaborn-v0_8-darkgrid")

plt.rcParams.update({
    "font.size": 12,
    "savefig.dpi": 150
})
```

For Jupyter Notebook
```python
matplotlib.style.use("dark_background")
matplotlib.rcParams.update({
    "axes.grid" : True,
    "grid.color" : "gray",
    "grid.linestyle" : "--"
})
```

Use svg in notebook
```python
import matplotlib_inline
matplotlib_inline.backend_inline.set_matplotlib_formats("svg")
```


