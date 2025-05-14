import numpy as np

#import tensorflow as tf
import tensorflow.compat.v1 as tf
tf.disable_v2_behavior()

from tensorflow.python.ops.array_grad import _TileGrad
from tensorflow.python.framework import ops


def shape(x):
    if isinstance(x, tf.Tensor):
        return x.get_shape().as_list()
    return np.shape(x)


@ops.RegisterGradient("TileDense")
def tile_grad_dense(op, grad):
    grad = tf.convert_to_tensor(grad)
    return _TileGrad(op, grad)
