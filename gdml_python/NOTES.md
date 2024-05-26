# About datasets

Datasets are given in npz or extended XYZ format.

In case of npz, it can be loaded using np.load.
It will return a descriptor object which cannot investigated directly.
It behaves like a loader or a dict. Some important keys are:
z, R, E, and F.

For example:

- dataset["E"][0] will return the energy of the first data

- dataset["R"][0] will return the coordinates of the first frame (as Natomsx3) array


# Windows vs Linux

On WSL2 and Windows, I got different results for alphas. The coefficients
are quite large numbers, the results are not matched perfectly both using
Cholesky and LU decomposition. The relative difference is about the the order 1e-6 or
1e-5 percent, so it is quite small if we look at the relative errors.
This difference might be attributed to rounding errors.
The K matrix and y coefficients are the same.

