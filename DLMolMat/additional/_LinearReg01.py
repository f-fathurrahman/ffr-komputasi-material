# %% [markdown]
# # Regresi linear sederhana
#
# Pada permasalahan *supervised learning*, kita biasanya berhadapan
# dengan masalah berikut. Diberikan $N$ pasangan data $\{(x_{i},t_{i})\}$ dengan
# $i=1,2,\ldots,N$ dan
#
# - $x$ adalah input atau masukan atau fitur, dan
# - $t$ adalah target atau output
#
# akan dicari output baru $t_{\mathrm{new}}$ dari input
# baru $x_{\mathrm{new}}$ yang diberikan.
#
# Model linear adalah salah satu model paling sederhana yang dapat digunakan.
# Pada model linear, diasumsikan $x$ dan $t$ memiliki hubungan linear.
# Secara matematis dapat dituliskan sebagai berikut
# $$
# t = f(x; w_0, w_1) = w_0 + w_1 x
# $$
# dengan $(w_0, w_1)$ adalah parameter dari model.

# %%
import numpy as np
