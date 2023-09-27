import matplotlib.pyplot as plt
plt.rcParams["savefig.dpi"] = 150

import rdkit.Chem

mol = rdkit.Chem.MolFromSmiles('C1CNCCC1C(=O)C')
img = rdkit.Chem.Draw.MolToImage(mol)

fig = plt.figure(figsize=(5,5))
ax = fig.gca()
ax.imshow(img)
ax.axis("off")
plt.savefig("IMG_test_mol.png")
