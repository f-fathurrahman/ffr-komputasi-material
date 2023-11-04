import matplotlib.pyplot as plt
import numpy as np

import matplotlib
matplotlib.style.use("seaborn-v0_8-darkgrid")

plt.rcParams.update({
    "savefig.dpi": 150
})


import pandas as pd
df = pd.read_csv("../../DATASET/iris.data.csv", header=None, encoding="utf-8")
# use of header=None is important
# Rename the columns
df.rename(columns={
    0 : "sepal-length", 1 :"sepal-width",
    2 : "petal-length", 3 :"petal-width",
    4 : "species"
}, inplace=True)

# Target class
class_all = df.iloc[:,4].values
# Get idx for each class
idx_set = class_all == "Iris-setosa"
idx_ver = class_all == "Iris-versicolor"
idx_vir = class_all == "Iris-virginica"
# Convert class to integer
class_idx = np.zeros(len(class_all), dtype=int)
class_idx[idx_set] = 0
class_idx[idx_ver] = 2
class_idx[idx_vir] = 1

# The features, choose only 2 features
# extract sepal length and petal length
X = df.iloc[:,[0, 2]].values

# plot data
plt.scatter(X[idx_set,0], X[idx_set,1], color="red", marker="o", label="Setosa")
plt.scatter(X[idx_ver,0], X[idx_ver,1], color="blue", marker="s", label="Versicolor")
plt.scatter(X[idx_vir,0], X[idx_vir,1], color="green", marker="^", label="Virginica")
plt.xlabel("Sepal length [cm]")
plt.ylabel("Petal length [cm]")
plt.legend(loc="upper left")
plt.savefig("IMG_iris_2.png")
plt.show()


plt.clf()
import seaborn as sns
sns.pairplot(df, hue="species")
plt.savefig("IMG_iris_all.png")
