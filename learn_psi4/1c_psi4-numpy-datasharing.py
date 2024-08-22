# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.16.2
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# %% [markdown]
# # Psi4 $\leftrightarrow$ NumPy Data Sharing
#
# The heart of the Psi4NumPy project its the ability to easily share and manipulate quantities in Python. While Psi4 offers the ability to manipulate most objects and perform tensor operations at the Python layer, it is often much easier to use the NumPy project, as its focus is on ease of use rather than optimal performance. Fortunately, Psi4 offers seemless integration with the NumPy framework. More details on the underlying functions can be found in the Psi4 [documentation](http://psicode.org/psi4manual/master/numpy.html).
#
# As before, let us start off with importing Psi4 and NumPy while also creating a random `5 x 5` NumPy array:

# %%
import psi4
import numpy as np

# Random number array
array = np.random.rand(5, 5)

# %% [markdown]
# Converting this to a Psi4 Matrix, which is an instance of the [`psi4.core.Matrix`](http://psicode.org/psi4manual/master/psi4api.html#psi4.core.Matrix 
# "Go to API") class, and back again is as simple as:

# %%
psi4_matrix = psi4.core.Matrix.from_array(array)
new_array = np.array(psi4_matrix)

print("Allclose new_array, array:", np.allclose(new_array, array))

# %% [markdown]
# ## Views
# Because both of these objects have the same in-memory data layout, the conversion is accomplished through the NumPy 
# [array_interface](https://docs.scipy.org/doc/numpy/reference/arrays.interface.html). This also opens the opportunity 
# to manipulate the Psi4 Matrix and Vector classes directly in memory.  To do this, we employ the `.np` attribute:

# %%
np.__version__

# %%
matrix = psi4.core.Matrix(3, 3)
print("Zero Psi4 Matrix:")
print(np.array(matrix))

matrix.np[:] = 1
print("\nMatrix updated to ones:")
print(np.array(matrix))

# %% [markdown]
# The `.np` attribute effectively returns a NumPy [view](http://scipy-cookbook.readthedocs.io/items/ViewsVsCopies.html). This view can then be manipulated as a conventional NumPy array and the underlying Psi4 Matrix data will be modified.
#
# <font color='red'>**Warning!** The following operation operation is incorrect and can potenitally lead to confusion:</font>

# %%
print(psi4.core.Matrix(3, 3).np)

# %% [markdown]
# While the above operation works about ~90% of the time, occasionally you will notice extremely large and small values. This is due to the fact that when you create the Psi4 Matrix and grab its view, the Psi4 Matrix is no longer bound to anything, and Python will attempt to "garbage collect" or remove the object. This sometimes happens *before* Python prints out the object so the NumPy view is pointing to a random piece of data in memory. A safe way to do this would be:

# %%
mat = psi4.core.Matrix(3, 3)
print(mat.np)

# or
print(np.asarray(psi4.core.Matrix(3, 3)))

# %% [markdown]
# Similar to the `.np` attribute, one can use `np.asarray` to create a NumPy view of a Psi4 object. Keep in mind that this is different than `np.array` which will copy the data.

# %%
mat = psi4.core.Matrix(3, 3)
mat_view = np.asarray(mat)

mat_view[:] = np.random.random(mat.shape)
print(mat.np)

# %% [markdown]
# Keep in mind that you must *update* this view using the `[]` syntax and not replace it (`=`). The following example should demonstrate the difference:

# %%
mat_view = np.zeros((3, 3))

# Mat is not updated as we replaced the mat_view with a new NumPy matrix.
print(mat.np)

# %% [markdown]
# ## Vector class
# Like the Psi4 Matrix class, the [`psi4.core.Vector`](http://psicode.org/psi4manual/master/psi4api.html#psi4.core.Vector "Go to API")
# class has similar accessors:

# %%
arr = np.random.rand(5)
vec = psi4.core.Vector.from_array(arr)
print(vec.np)
