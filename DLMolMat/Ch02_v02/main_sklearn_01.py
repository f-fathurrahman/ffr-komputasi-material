import numpy as np
from sklearn import linear_model
import matplotlib.pyplot as plt

GBL_a = 1.2
GBL_b = 3.1

def f_true(x1):
    return GBL_a*x1 + GBL_b

x = np.linspace(-5.0, 5.0, 10)
y = f_true(x)

# Cara manual
#Ndata = x.shape[0]
#X = np.zeros((Ndata,2))
#X[:,0] = 1.0
#X[:,1] = xs
#w = np.linalg.inv(X.T @ X) @ X.T @ y # bisa juga menggunakan SVD

# Inisialisasi model
model = linear_model.LinearRegression()
X = x.reshape(-1,1)
model.fit(X, y) # train
print("coef = ", model.coef_, " true param = ", GBL_a)
print("intercept = ", model.intercept_, " true param = ", GBL_b)



