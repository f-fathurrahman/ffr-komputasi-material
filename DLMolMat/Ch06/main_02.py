import pandas as pd

import tensorflow as tf
#tf.keras.backend.set_floatx("float64")

import numpy as np

import matplotlib.pyplot as plt
import matplotlib

matplotlib.style.use("seaborn-v0_8-darkgrid")
plt.rcParams.update({
    "font.size": 12,
    "savefig.dpi": 150
})


# Original: https://dataverse.harvard.edu/api/access/datafile/3407241
soldata = pd.read_csv("../DATASET/curated-solubility-dataset.csv")
features_start_at = list(soldata.columns).index("MolWt")
feature_names = soldata.columns[features_start_at:]
# standardize the features
soldata[feature_names] -= soldata[feature_names].mean()
soldata[feature_names] /= soldata[feature_names].std()

# Prepare data for Keras
full_data = tf.data.Dataset.from_tensor_slices(
    (soldata[feature_names].values, soldata["Solubility"].values)
)
N = len(soldata)
test_N = int(0.1 * N)
test_data = full_data.take(test_N).batch(32)
train_data = full_data.skip(test_N).batch(32)


# Build a neural network

# Hidden layers, only need to define output dimension
layer1 = tf.keras.layers.Dense(32, activation="tanh")
layer2 = tf.keras.layers.Dense(32, activation="tanh")

# Last layer. We want to output one number
output_layer = tf.keras.layers.Dense(1)

# Now put the layers into a sequential model
model = tf.keras.Sequential()
model.add(layer1)
model.add(layer2)
model.add(output_layer)

# Try our model on several datapoints
print("Try evaluation data:")
res = model(soldata[feature_names].values[:3])
print(res)
# Note that this is evaluated using current value of parameters in the model

# Prepare the model for training
model.compile(optimizer="Adam", loss="mean_squared_error")

# Train the model
model.fit(train_data, epochs=50)



def evaluate_performance(model, x, y, filename="IMG_01.png", title="Model Performance"):
    # get model predictions on test data and get labels
    # squeeze to remove extra dimensions
    yhat = np.squeeze(model.predict(x))

    plt.clf()
    plt.plot(y, yhat, ".")
    plt.plot(y, y, "-") # linear line
    plt.xlabel("Measured Solubility $y$")
    plt.ylabel("Predicted Solubility $\hat{y}$")
    plt.text(
        min(y) + 1,
        max(y) - 2,
        f"correlation = {np.corrcoef(y, yhat)[0,1]:.3f}",
    )
    plt.text(
        min(y) + 1,
        max(y) - 3,
        f"loss = {np.sqrt(np.mean((y - yhat)**2)):.3f}",
    )

    plt.gca().set_aspect("equal", "box")
    plt.title(title)
    plt.savefig(filename)
    plt.show()


test_y = soldata["Solubility"].values[:test_N]
evaluate_performance(model, test_data, test_y,
    filename="IMG_sol_02_test.png", title="Model Performance on Test Data")


train_y = soldata["Solubility"].values[test_N:]
evaluate_performance(model, train_data, train_y,
    filename="IMG_sol_02_train.png", title="Model Performance on Train Data")

