# %% [markdown]
# # Overfitting

# %% [markdown]
# Dataset biasanya dibagi menjadi setidaknya dua subset:
# - data latih (train dataset)
# - data uji (test dataset)

# %% [markdown]
# Pada deep learning biasanya juga digunakan satu dataset lagi:
# - data validasi (validation dataset)

# %% [markdown]
# Data latih digunakan untuk mendapatkan parameter model, sedangkan data uji digunakan untuk mengevaluasi performa model (apakah model bekerja dengan baik atau tidak).

# %% [markdown]
# Banyak framework seperti TensorFlow, PyTorch, dan Scikit Learn memiliki tools sendiri untuk menangani dataset.

# %%
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import jax.numpy as jnp
import jax

# %%
soldata = pd.read_csv("../DATASET/curated-solubility-dataset.csv")
features_start_at = list(soldata.columns).index("MolWt")
feature_names = soldata.columns[features_start_at:]


# %%
len(soldata)

# %% [markdown]
# Hanya 50 data yang digunakan

# %%
# Get 50 points and split into train/test
sample = soldata.sample(50, replace=False)
train = sample[:25] # data latih
test = sample[25:] # data uji

# %% [markdown]
# Standardisasi fitur (secara manual):

# %%
test[feature_names] -= train[feature_names].mean()
test[feature_names] /= train[feature_names].std()
train[feature_names] -= train[feature_names].mean()
train[feature_names] /= train[feature_names].std()

# convert from pandas dataframe to numpy arrays
x = train[feature_names].values
y = train["Solubility"].values
test_x = test[feature_names].values
test_y = test["Solubility"].values

# %%
x.mean()

# %%
x.std(ddof=1)


# %%
# define our loss function
def loss(w, b, x, y):
    return jnp.mean( (y - jnp.dot(x, w) - b)**2)

loss_grad = jax.grad(loss, (0, 1))
# Inisialisasi parameter secara acak
w = np.random.normal(size=x.shape[1])
b = 0.0
loss_grad(w, b, x, y)

# %%
loss_progress = []
test_loss_progress = []
eta = 0.05
# belum menggunakan batching
for i in range(2000):
    grad = loss_grad(w, b, x, y)
    w -= eta * grad[0]
    b -= eta * grad[1]
    loss_progress.append(loss(w, b, x, y))
    test_loss_progress.append(loss(w, b, test_x, test_y))

# %%
plt.plot(loss_progress, label="Training Loss")
plt.plot(test_loss_progress, label="Testing Loss")
plt.legend()

# %% [markdown]
# Salah satu yang bisa dilakukan adalah: "early stopping", sebelum test lost naik, training dihentikan.

# %%
yhat = x @ w + b
plt.plot(y, y, ":", linewidth=0.2)
plt.plot(y, x @ w + b, "o")
plt.xlim(min(y), max(y))
plt.ylim(min(y), max(y))
plt.text(min(y) + 1, max(y) - 2, f"correlation = {np.corrcoef(y, yhat)[0,1]:.3f}")
plt.text(min(y) + 1, max(y) - 3, f"loss = {np.sqrt(np.mean((y - yhat)**2)):.3f}")
plt.title("Training Data")

# %%
yhat = test_x @ w + b
plt.plot(test_y, test_y, ":", linewidth=0.2)
plt.plot(test_y, yhat, "o")
plt.xlim(min(test_y), max(test_y))
plt.ylim(min(test_y), max(test_y))
plt.text(
    min(test_y) + 1,
    max(test_y) - 2,
    f"correlation = {np.corrcoef(test_y, yhat)[0,1]:.3f}",
)
plt.text(
    min(test_y) + 1,
    max(test_y) - 3,
    f"loss = {np.sqrt(np.mean((test_y - yhat)**2)):.3f}",
)
plt.title("Testing Data")

# %% [markdown]
# # Data sintetik

# %%
# generate data from polynomial
N = 20
syn_x = np.linspace(-3, 3, N)
# create feature matrix
syn_features = np.vstack([syn_x**3, syn_x**2, syn_x, np.ones_like(syn_x)]).T
syn_labels = syn_x**3 - syn_x**2 + syn_x - 1 # true model

# %%
# split data into train/test
indices = list(range(0, N // 4)) + list(range(3 * N // 4, N))
test_indices = list(range(N // 4, 3 * N // 4))

# %%
indices

# %%
test_indices

# %%
train_x = syn_features[indices]
train_y = syn_labels[indices]
test_x = syn_features[test_indices]
test_y = syn_labels[test_indices]

# %%
# fit using numpy least squares method.
w, *_ = np.linalg.lstsq(train_x, train_y, rcond=-1)

# %%
w

# %%
# plotting code
plt.plot(syn_x[indices], train_y, "o", label="Training labels")
plt.plot(syn_x[test_indices], test_y, "o", label="Testing labels")
plt.ylim(-40, 40)
plt.plot(syn_x, jnp.dot(syn_features, w), label="Fit Model")
plt.plot(syn_x, syn_labels, "--", label="Ground Truth")
plt.text(0, -20, f"Training Loss {loss(w,0,train_x, train_y):.2f}")
plt.text(0, -30, f"Testing Loss {loss(w,0, test_x, test_y):.2f}")
plt.legend()
plt.title("No Noise, Perfect Features")
plt.show()

# %% [markdown]
# Tambah noise:

# %%
train_y = train_y + 0.01*np.random.normal(scale=1.0, size=train_y.shape)
test_y = test_y + 0.01*np.random.normal(scale=1.0, size=test_y.shape)

# %%
w, *_ = np.linalg.lstsq(train_x, train_y, rcond=-1)
plt.plot(syn_x[indices], train_y, "o", label="Training labels")
plt.plot(syn_x[test_indices], test_y, "o", label="Testing labels")
plt.ylim(-40, 40)
plt.plot(syn_x, jnp.dot(syn_features, w), label="Fit Model")
plt.plot(syn_x, syn_labels, "--", label="Ground Truth")
plt.text(0, -20, f"Training Loss {loss(w,0,train_x, train_y):.2f}")
plt.text(0, -30, f"Testing Loss {loss(w,0, test_x, test_y):.2f}")
plt.legend()
plt.title("Noise, Perfect Features")

# %%
syn_features = np.vstack([syn_x**i for i in range(7)]).T

# %%
train_x = syn_features[indices]
test_x = syn_features[test_indices]
test_y = syn_labels[test_indices]

w, *_ = np.linalg.lstsq(train_x, train_y, rcond=-1)
plt.plot(syn_x[indices], train_y, "o", label="Training labels")
plt.plot(syn_x[test_indices], test_y, "o", label="Testing labels")
plt.ylim(-40, 40)
plt.plot(syn_x, jnp.dot(syn_features, w), label="Fit Model")
plt.plot(syn_x, syn_labels, "--", label="Ground Truth")
plt.text(0, -20, f"Training Loss {loss(w,0,train_x, train_y):.2f}")
plt.text(0, -30, f"Testing Loss {loss(w,0, test_x, test_y):.2f}")
plt.legend(loc="upper left")
plt.title("Noise, Extra Features")

# %%
