# %% [markdown]
# ## Operasi array

# %%
import tensorflow as tf



# %% [markdown]
# Bilangan acak

# %%
A = tf.random.uniform([3, 4])

# %%
B = tf.random.normal([4, 3])

# %%
type(A)

# %% [markdown]
# Menggunakan `tf.Variable`:

# %%
a = tf.Variable([1.0, 2.0, 3.0, 4.0])

# %%
type(a)

# %%
a.get_shape()

# %%
A * a

# %%
a * A

# %% [markdown]
# Variabel lain

# %%
b = tf.Variable([[1.0], [2.0], [3.0]])

# %%
b.get_shape()

# %%
A * b

# %%
Bb = tf.linalg.matmul(B, b)
Bb

# %%
type(Bb)

# %%
Br = tf.reshape(B, (3,4))

# %%
Br, B

# %% [markdown]
# Matrix multiplication using `tf.Variable`

# %%
C = tf.Variable([[1.0, 2.0], [3.0, 4.0]])
x = tf.Variable([[1.0], [2.0]])

# %%
Cx = tf.linalg.matmul(C, x)
Cx, type(Cx)

# %% [markdown]
# Transpose:

# %%
Bt = tf.transpose(B)

# %%
Bt, B

# %%
tf.linalg.matmul(Bt, B)

# %%
Bt.shape

# %%
B.shape

# %% [markdown]
# Bilangan kompleks

# %%
Z1 = tf.Variable([
    [1 + 1j, 1 + 2.6j],
    [1 - 2.5j, 3 + 2.3j]
])

# %%
Z1.dtype

# %%
Z1.shape

# %%
tf.transpose(Z1)

# %%
Z1, tf.math.conj(Z1)

# %%
tf.transpose( tf.math.conj(Z1) )

# %%
