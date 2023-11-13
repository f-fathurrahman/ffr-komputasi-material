import numpy as np
import tensorflow as tf

x = np.random.normal(size=[1,4]).astype("float32")

init = tf.keras.initializers.RandomNormal(seed=1234)

w1 = tf.Variable(init(shape=[4,3]))
b1 = tf.Variable(init(shape=[1,3]))

w2 = tf.Variable(init(shape=[3,2]))
b2 = tf.Variable(init(shape=[1,2]))

@tf.function
def forward(x, W, b, activation):
    return activation(tf.matmul(x, W) + b)

h = forward(x, w1, b1, tf.nn.sigmoid)
y = forward(h, w2, b2, tf.nn.softmax)

print(y)


