import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

# RDkit stuffs are not included

# Original: https://dataverse.harvard.edu/api/access/datafile/3407241
soldata = pd.read_csv("../DATASET/curated-solubility-dataset.csv")

# Plot histogram
sns.distplot(soldata.Solubility)
plt.show()


# Plot again
features_start_at = list(soldata.columns).index("MolWt")
feature_names = soldata.columns[features_start_at:]

fig, axs = plt.subplots(nrows=5, ncols=4, sharey=True, figsize=(12, 8), dpi=300)
axs = axs.flatten()  # so we don't have to slice by row and column
for i, n in enumerate(feature_names):
    ax = axs[i]
    ax.scatter(
        soldata[n], soldata.Solubility, s=6, alpha=0.4, color=f"C{i}"
    )  # add some color
    if i % 4 == 0:
        ax.set_ylabel("Solubility")
    ax.set_xlabel(n)

# hide empty subplots
for i in range(len(feature_names), len(axs)):
    fig.delaxes(axs[i])

plt.tight_layout()
plt.savefig("IMG_grid_plot_01.png", dpi=200)
plt.savefig("IMG_grid_plot_01.pdf")
#plt.show()
