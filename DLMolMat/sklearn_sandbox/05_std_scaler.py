import numpy as np
from sklearn.preprocessing import PolynomialFeatures, StandardScaler

X = np.linspace(-1.0, 1.0, 5)[:,np.newaxis]

polyN = PolynomialFeatures(10)
std_scaler = StandardScaler()

polyN.fit(X)
X1 = polyN.transform(X)

std_scaler.fit(X1)
X1 = std_scaler.transform(X1)
print("X1 = ")
print(X1)

# different X
Xnew = np.linspace(-1.0, 1.0, 6)[:,np.newaxis]
# use the same transform, skip fitting step
X2 = polyN.transform(Xnew)
X2 = std_scaler.transform(X2)

print("X2 = ")
print(X2)
