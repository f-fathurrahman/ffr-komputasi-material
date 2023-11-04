import matplotlib.pyplot as plt
import numpy as np

import matplotlib
matplotlib.style.use("seaborn-v0_8-darkgrid")

plt.rcParams.update({
    "savefig.dpi": 150
})


import pandas as pd
df = pd.read_csv("../../DATASET/iris.data.csv", header=None, encoding='utf-8')
# use of header=None is important

# select setosa and versicolor (first 100 rows)
y = df.iloc[0:100, 4].values
y = np.where(y == 'Iris-setosa', 0, 1)
# Iris-setosa will be given 0 and Iris-versicolor will be given 1

# The features, choose only 2 features
# extract sepal length and petal length
X = df.iloc[0:100, [0, 2]].values

# plot data
plt.scatter(X[:50, 0], X[:50, 1], color='red', marker='o', label='Setosa')
plt.scatter(X[50:100, 0], X[50:100, 1], color='blue', marker='s', label='Versicolor')
plt.xlabel('Sepal length [cm]')
plt.ylabel('Petal length [cm]')
plt.legend(loc='upper left')
plt.savefig("IMG_iris_setosa_versicolor.png")
plt.show()

