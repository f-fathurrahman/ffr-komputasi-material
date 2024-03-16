import pandas as pd
import numpy as np
from sklearn import linear_model
import matplotlib.pyplot as plt

soldata = pd.read_csv("../DATASET/curated-solubility-dataset.csv")
features_start_at = list(soldata.columns).index("MolWt")
feature_names = soldata.columns[features_start_at:]

# convert data into features, labels
features = soldata.loc[:, feature_names].values
labels = soldata.Solubility.values

model = linear_model.LinearRegression() # bandingkan dengan MLPRegressor

model.fit(features, labels)
predicted_labels = model.predict(features)

# Buat plot paritas
plt.scatter(labels, predicted_labels, s=4, alpha=0.7)
plt.xlabel("Measured Solubility $y$")
plt.ylabel("Predicted Solubility $\hat{y}$")
plt.xlim(min(labels), max(labels))
plt.ylim(min(labels), max(labels))
plt.gca().set_aspect("equal", "box")
plt.savefig("IMG_parity_plot.pdf")
plt.show()
