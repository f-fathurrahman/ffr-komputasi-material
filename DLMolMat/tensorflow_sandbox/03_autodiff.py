# %%
import tensorflow as tf

# %% [markdown]
# # Diferensiasi Otomatis (Automatic Diferentiation)

# %% [markdown]
# Untuk menghitungan turunan fungsi secara otomatis, TensorFlow perlu merekam proses perhitungan apasaja yang dilakukan.

# %% [markdown]
# TensorFlow provides the tf.GradientTape API for automatic differentiation; that is, computing the gradient of a computation with respect to some inputs, usually tf.Variables. TensorFlow "records" relevant operations executed inside the context of a tf.GradientTape onto a "tape". TensorFlow then uses that tape to compute the gradients of a "recorded" computation using reverse mode differentiation.

# %%
x = tf.Variable(3.0)

# %% [markdown]
# Lakukan perhitungan, rekam dengan menggunakan `tf.GradientTape`

# %%
with tf.GradientTape() as tape:
  y = 0.4*x**2


# %%
y

# %%
dy_dx = tape.gradient(y, x)

# %%
dy_dx

# %% [markdown]
# Fungsi atau perhitungan yang dilakukan adalah $y = 0.4x^2$.
#
# Turunannya adalah: $y' = 0.8x$.
#
# Dievaluasi pada $x = 3$ memberikan:

# %%
0.4*3**2

# %%
0.8*3

# %% [markdown]
# # Pada vektor/tensor

# %%
w = tf.Variable(tf.random.normal((3, 2)), name='w')
b = tf.Variable(tf.zeros(2, dtype=tf.float32), name='b')
x = [[1., 2., 3.]]

with tf.GradientTape(persistent=True) as tape:
  y = x @ w + b
  loss = tf.reduce_mean(y**2)

# %%
[dl_dw, dl_db] = tape.gradient(loss, [w, b])

# %%
dl_dw

# %%
dl_db

# %%
