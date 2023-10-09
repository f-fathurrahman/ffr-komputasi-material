import numpy as np
from scipy.signal import convolve2d

import matplotlib.pyplot as plt
from skimage.io import imread, imsave
from skimage.util import img_as_float32, img_as_ubyte

img = imread("../DATASET/The_Cat.jpg")

plt.figure(figsize=(6,10))
plt.imshow(img)
plt.show()


