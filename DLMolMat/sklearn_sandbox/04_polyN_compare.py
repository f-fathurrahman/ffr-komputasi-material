import numpy as np
from sklearn.preprocessing import PolynomialFeatures

X = np.linspace(-1.0, 1.0, 5)[:,np.newaxis]

polyN = PolynomialFeatures(3)
polyN.fit(X)
X1 = polyN.transform(X)
print("X1 = ")
print(X1)

# different X
Xnew = np.linspace(-1.0, 1.0, 6)[:,np.newaxis]
# use the same transform, skip fitting step
X2 = polyN.transform(Xnew)
print("X2 = ")
print(X2)
