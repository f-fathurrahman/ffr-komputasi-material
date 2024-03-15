import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import rdkit, rdkit.Chem

soldata = pd.read_csv("../DATASET/curated-solubility-dataset.csv")
soldata.head()

#sns.distplot(soldata.Solubility) # deprecated
sns.histplot(soldata.Solubility)
plt.savefig("IMG_histplot.png", dpi=200)
plt.savefig("IMG_histplot.pdf")
plt.show()


# get 3 lowest and 3 highest solubilities
soldata_sorted = soldata.sort_values("Solubility")
extremes = pd.concat([soldata_sorted[:3], soldata_sorted[-3:]])

# We need to have a list of strings for legends
legend_text = [
    f"{x.ID}: solubility = {x.Solubility:.2f}" for x in extremes.itertuples()
]

# now plot them on a grid
extreme_mols = [rdkit.Chem.MolFromInchi(inchi) for inchi in extremes.InChI]
images = rdkit.Chem.Draw.MolsToGridImage(
    extreme_mols, molsPerRow=3, subImgSize=(500, 500), legends=legend_text
)
images.save("IMG_extremes.png")



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
plt.savefig("IMG_feature_correlation.png")
plt.show()

