import numpy as np

batch_size = 4
N = 16  # number of data points
# compute how much data fits nicely into a batch and drop extra data
new_N = N // batch_size * batch_size

# Dummy data
feature_dim = 3
features = np.random.rand(N,feature_dim)
labels = np.random.rand(N)

# the -1 means that numpy will compute
# what that dimension should be
batched_features = features[:new_N].reshape((-1, batch_size, feature_dim))
batched_labels = labels[:new_N].reshape((-1, batch_size))
indices = np.arange(new_N // batch_size)

print("indices = ", indices)