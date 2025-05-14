#import tensorflow as tf
import tensorflow.compat.v1 as tf
tf.disable_v2_behavior()


def _softplus(x):
    return tf.log1p(tf.exp(x))


def shifted_softplus(x):
    """
    Softplus nonlinearity shifted by -log(2) such that shifted_softplus(0.) = 0.

    y = log(0.5e^x + 0.5)

    """
    y = tf.where(x < 14., _softplus(tf.where(x < 14., x, tf.zeros_like(x))), x)
    return y - tf.log(2.)
