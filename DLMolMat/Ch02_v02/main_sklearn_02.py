import numpy as np
from sklearn import linear_model
import matplotlib.pyplot as plt

GBL_a1 = 1.2
GBL_a2 = 9.9
GBL_b = 3.1

def f_true(x1, x2):
    return GBL_a1*x1 + GBL_a2*x2 + GBL_b

Ndata = 10
#x1 = np.random.rand(Ndata)
#x2 = np.random.rand(Ndata)
x1 = np.linspace(-1.0, 1.0, Ndata)
x2 = np.linspace(-2.1, 2.0, Ndata) #+ np.random.randn(Ndata)
y = f_true(x1, x2) #+ np.random.randn()

# Cara manual
X = np.zeros((Ndata,3))
X[:,0] = 1.0
X[:,1] = x1
X[:,2] = x2
GG = X.T @ X
print("Condition number of X: ", np.linalg.cond(X))
print("Condition number of Gramian: ", np.linalg.cond(GG))
w = np.linalg.inv(GG) @ X.T @ y # bisa juga menggunakan SVD
print(w)


print("Using sklearn LinearRegression: ")
# Inisialisasi model
model = linear_model.LinearRegression()
X = np.zeros((Ndata,2))
X[:,0] = x1
X[:,1] = x2
model.fit(X, y)
print("coef = ", model.coef_, " true param = ", [GBL_a1, GBL_a2])
print("intercept = ", model.intercept_, " true param = ", GBL_b)
