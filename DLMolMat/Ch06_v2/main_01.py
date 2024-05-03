# %% [markdown]
# # Chapter 6

# %%
import pandas as pd
import matplotlib.pyplot as plt
import tensorflow as tf
import numpy as np

# %%
tf.config.list_physical_devices("CPU")


# %%
soldata = pd.read_csv("../DATASET//curated-solubility-dataset.csv")
features_start_at = list(soldata.columns).index("MolWt")
feature_names = soldata.columns[features_start_at:]
# standardize the features
soldata[feature_names] -= soldata[feature_names].mean()
soldata[feature_names] /= soldata[feature_names].std()

# %%
len(feature_names)

# %% [markdown]
# - Convert from numpy array to tf.Dataset
# - Split data into training and test dataset

# %%
full_data = tf.data.Dataset.from_tensor_slices(
    (soldata[feature_names].values, soldata["Solubility"].values)
)
N = len(soldata)
test_N = int(0.1 * N) # 10% of all data as test data
test_data = full_data.take(test_N).batch(16)
train_data = full_data.skip(test_N).batch(16) # 16: batch size

# %%
test_N

# %%
N

# %% [markdown]
# Build the network

# %%
# our hidden layer
# We only need to define the output dimension 
hidden_layer1 = tf.keras.layers.Dense(32, activation="tanh")
hidden_layer2 = tf.keras.layers.Dense(16, activation="tanh")

# %%
# Last layer - which we want to output one number
# the predicted solubility.
output_layer = tf.keras.layers.Dense(1)

# %%
# Now we put the layers into a sequential model
model = tf.keras.Sequential()
# input layer secara default dapat dihandle oleh Keras
model.add(hidden_layer1)
model.add(hidden_layer2)
model.add(output_layer)

# %%
xx = soldata[feature_names].values[:1]
xx.shape

# %%
# Try out our model on first few datapoints
model(soldata[feature_names].values[:1])

# %%
model.summary()

# %%
17*32 + 32

# %%
17*50 + 50 # 17: feature number, 50 hidden node, 50 bias

# %%
model.compile(optimizer="SGD", loss="mean_squared_error")

# %%
model.fit(train_data, epochs=50)

# %%
tf.keras.layers.

# %%
yhat = np.squeeze(model.predict(test_data))
yhat.shape

# %%
# get model predictions on test data and get labels
# squeeze to remove extra dimensions
test_y = soldata["Solubility"].values[:test_N]

# %%
np.sqrt(np.mean((test_y - yhat)**2)) # root mean squared error

# %%
import matplotlib
matplotlib.style.use("dark_background")
matplotlib.rcParams.update({
    "axes.grid" : True,
    "grid.color" : "gray",
    "grid.linestyle" : "--"
})

# %%
plt.plot(test_y, yhat, ".")
plt.plot(test_y, test_y, "-")
plt.xlabel("Measured Solubility $y$")
plt.ylabel("Predicted Solubility $\hat{y}$")
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
plt.show()

# %%
plt.plot(test_y, yhat, ".")
plt.plot(test_y, test_y, "-")
plt.xlabel("Measured Solubility $y$")
plt.ylabel("Predicted Solubility $\hat{y}$")
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
plt.show()

# %%
