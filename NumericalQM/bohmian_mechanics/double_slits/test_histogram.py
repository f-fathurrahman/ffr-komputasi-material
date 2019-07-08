#!/usr/bin/env python

import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import scipy.stats

# example data
mu = 0.0 # mean of distribution
sigma = 1.0 # standard deviation of distribution
x = mu + sigma * np.random.randn(10000)

num_bins = 20
# the histogram of the data
n, bins, patches = plt.hist(x, num_bins, density=True, facecolor="blue", alpha=0.5)

# add a "best fit" line
#y = mlab.normpdf(bins, mu, sigma)
y = scipy.stats.norm.pdf(bins, mu, sigma)
plt.plot(bins, y, "r--")
plt.xlabel("Variable")
plt.ylabel("Probability")
plt.title(r"Normal distribution: $\mu=0$, $\sigma=1$")

# Tweak spacing to prevent clipping of ylabel
plt.subplots_adjust(left=0.15)
plt.savefig("TEMP_histogram.pdf")