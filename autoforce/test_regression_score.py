import sys
from ase.io import read
import numpy as np

from regression_scores import get_energy_and_forces, cd, maxe, mae, rmse

try:
    r = sys.argv[3]
except IndexError:
    r = '::'
data = read(sys.argv[1], r)
targets = read(sys.argv[2], r)

d_e, d_f, d_n = get_energy_and_forces(data)
t_e, t_f, t_n = get_energy_and_forces(targets)

assert (d_n == t_n).all()
quant = len(data)
assert len(targets) == quant

def _print(d, t):
    s_d = np.float64(np.sqrt(d.var()))
    s_t = np.float64(np.sqrt(t.var()))
    print(f'\tstdev:  {s_d} (p), {s_t} (t)')
    print(f'\tmaxe:   {maxe(d, t)}')
    print(f'\tmae:    {mae(d, t)}')
    print(f'\trmse:   {rmse(d, t)}')
    print(f'\tcd:     {cd(d, t)}')

print(f'predictions:       {sys.argv[1]}')
print(f'targets:           {sys.argv[2]}')
print(f'number of samples: {quant}')

for n, d, t in [('energy', d_e, t_e), ('forces', d_f, t_f)]:
    print(f'\n{n}:')
    _print(d, t)
    if n == 'forces':
        zset = np.unique(d_n)
        if len(zset) > 1:
            for z in zset:
                i = d_n == z
                s = i.sum()//3
                print(f'\nforces on (atomic number): {z}  (size= 3*{s})')
                _print(d[i], t[i])
