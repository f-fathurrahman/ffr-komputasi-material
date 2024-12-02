Descriptor files are quite big: for example SiO2 `~28GB`


```text
<redacted>pyxtal_ff\models_optimizers_lbfgs.py:264: UserWarning: This overload of add_ is deprecated:
        add_(Number alpha, Tensor other)
Consider using one of the following signatures instead:
        add_(Tensor other, *, Number alpha) (Triggered internally at ..\torch\csrc\utils\python_arg_parser.cpp:1578.)
  p.data.add_(step_size, update[offset:offset + numel].view_as(p.data))
```

