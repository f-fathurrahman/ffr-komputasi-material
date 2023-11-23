# %%
import tensorflow as tf

# %% [markdown]
# # Diferensiasi Otomatis (Automatic Diferentiation)

# %% [markdown]
# ## Intro

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
# ## Pada vektor/tensor

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

# %% [markdown]
# ## Gradien dari suatu model

# %%
layer = tf.keras.layers.Dense(2, activation='relu')
x = tf.Variable([[1., 2., 3.]])

# %%
with tf.GradientTape() as tape:
  # Forward pass
  y = layer(x)
  loss = tf.reduce_mean(y**2)

# Calculate gradients with respect to every trainable variable
grad = tape.gradient(loss, layer.trainable_variables)

# %%
grad

# %%
for var, g in zip(layer.trainable_variables, grad):
  print(f'{var.name}, shape: {g.shape}')


# %% [markdown]
# ## Mengontrol perilaku GradientTape

# %%
# A trainable variable
x0 = tf.Variable(0.5, name='x0')
# Not trainable
x1 = tf.Variable(3.0, name='x1', trainable=True)
# Not a Variable: A variable + tensor returns a tensor.
x2 = tf.Variable(2.0, name='x2') + 1.0
# Not a variable
x3 = tf.constant(3.0, name='x3')

with tf.GradientTape() as tape:
  y = (x0**3) + (x1**3) + (x2**3)

grad = tape.gradient(y, [x0, x1, x2, x3])

for g in grad:
  print(g)

# %% [markdown]
# You can list the variables being watched by the tape using the GradientTape.watched_variables method:

# %%
[var.name for var in tape.watched_variables()]

# %% [markdown]
# To record gradients with respect to a tf.Tensor, you need to call GradientTape.watch(x):

# %%
x = tf.constant(3.0)
with tf.GradientTape() as tape:
  tape.watch(x)
  y = x**2

# dy = 2x * dx
dy_dx = tape.gradient(y, x)
print(dy_dx.numpy())

# %%
[var.name for var in tape.watched_variables()]

# %% [markdown]
# Conversely, to disable the default behavior of watching all tf.Variables, set watch_accessed_variables=False when creating the gradient tape. This calculation uses two variables, but only connects the gradient for one of the variables:

# %%
x0 = tf.Variable(2.0)
x1 = tf.Variable(10.0)

with tf.GradientTape(watch_accessed_variables=False) as tape:
    tape.watch(x1)
    y0 = tf.math.sin(x0)
    y1 = tf.nn.softplus(x1)
    y = y0 + y1
    ys = tf.reduce_sum(y)

# dys/dx1 = exp(x1) / (1 + exp(x1)) = sigmoid(x1)
grad = tape.gradient(ys, {'x0': x0, 'x1': x1})

grad

# %%
