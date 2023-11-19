# %% [markdown]
# Eager execution


# %%
import platform
import tensorflow as tf

# %%
platform.python_version()

# %%
tf.__version__

# %%
tf.executing_eagerly()


# %%
def eager_func(a, b):
    return a + b


# %%
@tf.function
def graph_func(a, b):
    return a + b


# %%
a = tf.random.uniform([5,5])
b = tf.random.uniform([5,5])

# %%
# %timeit -n 100 eager_func(a, b)

# %%
# %timeit -n 100 graph_func(a, b)

# %% [markdown]
# More complex calculations

# %%
tf.keras.backend.clear_session()

# %%
from tensorflow.keras import Input, Model, Sequential
from tensorflow.keras.layers import Flatten, Dense

# %%
model = Sequential([
    Input(shape=(28,28)),
    Flatten(),
    Dense(256, "relu"),
    Dense(128, "relu"),
    Dense(256, "relu"),
    Dense(10, "softmax")
])

# %%
model.summary()

# %%
# Dummy input with MNIST image size
X = tf.random.uniform([1000, 28, 28])

# %%
type(X), X.dtype

# %%
import timeit

# %%
timeit.timeit(lambda: model(X, training=False), number=100)

# %%
timeit.timeit(lambda: model(X, training=False), number=100)

# %% [markdown]
# Using TensorFlow graph, compiled function by `tf.function`:

# %%
graph_model = tf.function(model)

# %%
timeit.timeit(lambda: graph_model(X, training=False), number=100)

# %%
timeit.timeit(lambda: graph_model(X, training=False), number=100)

# %%
