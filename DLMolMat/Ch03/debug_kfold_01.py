import numpy as np

k = 10 # number of folds
N = 110
# make indices for the k segments
splits = list(range(0, N + N//k, N//k))
print("splits = ", splits)
error = []
# slice out segments
for i in range(k):
    print("start splits[i] = ", splits[i], " splits[i+1] = ", splits[i+1])
