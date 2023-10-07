# Original
# https://github.com/che-shr-cat/JAX-in-Action/blob/main/Chapter-2/JAX_in_Action_Chapter_2_MNIST_MLP_Pure_JAX.ipynb

import tensorflow as tf
import tensorflow_datasets as tfds

data_dir = '../DATASET/tfds'
# as_supervised=True gives us the (image, label) as a tuple instead of a dict
data, info = tfds.load(
    name="mnist",
    data_dir=data_dir,
    as_supervised=True,
    with_info=True,
    download=False
)

data_train = data['train']
data_test  = data['test']

# XXX: How to work with tfds Dataset?



import numpy as np
import matplotlib.pyplot as plt
plt.rcParams["figure.figsize"] = [10, 5]

ROWS = 3
COLS = 10
i = 0
fig, ax = plt.subplots(ROWS, COLS)
for image, label in data_train.take(ROWS*COLS):
    ax[int(i/COLS), i%COLS].axis('off')
    ax[int(i/COLS), i%COLS].set_title(str(label.numpy()))
    ax[int(i/COLS), i%COLS].imshow(np.reshape(image, (28,28)), cmap='gray')
    i += 1

plt.show()




HEIGHT = 28
WIDTH  = 28
CHANNELS = 1
NUM_PIXELS = HEIGHT * WIDTH * CHANNELS
NUM_LABELS = info.features['label'].num_classes


import jax.numpy as jnp

def preprocess(img, label):
    """Resize and preprocess images."""
    return (tf.cast(img, tf.float32)/255.0), label



data_train_vis = data_train.map(preprocess)

i = 0
fig, ax = plt.subplots(ROWS, COLS)
for image, label in data_train_vis.take(ROWS*COLS):
    ax[int(i/COLS), i%COLS].axis('off')
    ax[int(i/COLS), i%COLS].set_title(str(label.numpy()))
    ax[int(i/COLS), i%COLS].imshow(np.reshape(image, (28,28)), cmap='gray')
    i += 1

plt.show()



train_data = tfds.as_numpy(data_train.map(preprocess).batch(32).prefetch(1))
test_data  = tfds.as_numpy(data_test.map(preprocess).batch(32).prefetch(1))


LAYER_SIZES = [28*28, 512, 10]
PARAM_SCALE = 0.01

import jax
import jax.numpy as jnp
from jax import grad, jit, vmap, value_and_grad
from jax import random
from jax.nn import swish, logsumexp, one_hot



def init_network_params(sizes, key=random.PRNGKey(0), scale=1e-2):
    """Initialize all layers for a fully-connected neural network with given sizes"""

    def random_layer_params(m, n, key, scale=1e-2):
        """A helper function to randomly initialize weights and biases of a dense layer"""
        w_key, b_key = random.split(key)
        return scale * random.normal(w_key, (n, m)), scale * random.normal(b_key, (n,))

    keys = random.split(key, len(sizes))
    return [random_layer_params(m, n, k, scale) for m, n, k in zip(sizes[:-1], sizes[1:], keys)]



params = init_network_params(LAYER_SIZES, random.PRNGKey(0), scale=PARAM_SCALE)




def predict(params, image):
    """Function for per-example predictions."""
    activations = image
    for w, b in params[:-1]:
        outputs = jnp.dot(w, activations) + b
        activations = swish(outputs)
    #
    final_w, final_b = params[-1]
    logits = jnp.dot(final_w, activations) + final_b
    return logits


random_flattened_image = random.normal(random.PRNGKey(1), (28*28*1,))
preds = predict(params, random_flattened_image)
print(preds.shape)



random_flattened_images = random.normal(random.PRNGKey(1), (32, 28*28*1))
try:
    preds = predict(params, random_flattened_images)
except TypeError as e:
    print(e)


batched_predict = vmap(predict, in_axes=(None, 0))

# `batched_predict` has the same call signature as `predict`
batched_preds = batched_predict(params, random_flattened_images)
print(batched_preds.shape)



jax.devices()


INIT_LR = 1.0
DECAY_RATE = 0.95
DECAY_STEPS = 5

def loss(params, images, targets):
    """Categorical cross entropy loss function."""
    logits = batched_predict(params, images)
    log_preds = logits - logsumexp(logits)
    # logsumexp trick https://gregorygundersen.com/blog/2020/02/09/log-sum-exp/
    return -jnp.mean(targets*log_preds)

@jit
def update(params, x, y, epoch_number):
    loss_value, grads = value_and_grad(loss)(params, x, y)
    lr = INIT_LR * DECAY_RATE ** (epoch_number / DECAY_STEPS)
    return [(w - lr * dw, b - lr * db)
          for (w, b), (dw, db) in zip(params, grads)], loss_value



num_epochs = 25

@jit
def batch_accuracy(params, images, targets):
  images = jnp.reshape(images, (len(images), NUM_PIXELS))
  predicted_class = jnp.argmax(batched_predict(params, images), axis=1)
  return jnp.mean(predicted_class == targets)

def accuracy(params, data):
  accs = []
  for images, targets in data:
    accs.append(batch_accuracy(params, images, targets))
  return jnp.mean(jnp.array(accs))

import time

for epoch in range(num_epochs):
    start_time = time.time()
    losses = []
    for x, y in train_data:
        x = jnp.reshape(x, (len(x), NUM_PIXELS))
        y = one_hot(y, NUM_LABELS)
        params, loss_value = update(params, x, y, epoch)
        losses.append(loss_value)
    epoch_time = time.time() - start_time

    start_time = time.time()
    train_acc = accuracy(params, train_data)
    test_acc = accuracy(params, test_data)
    eval_time = time.time() - start_time
    print("Epoch {} in {:0.2f} sec".format(epoch, epoch_time))
    print("Eval in {:0.2f} sec".format(eval_time))
    print("Training set loss {}".format(jnp.mean(jnp.array(losses))))
    print("Training set accuracy {}".format(train_acc))
    print("Test set accuracy {}".format(test_acc))



"""
# save the model
# https://github.com/google/jax/issues/2116

import pickle

model_weights_file = 'mlp_weights.pickle'

with open(model_weights_file, 'wb') as file:
    pickle.dump(params, file)


with open(model_weights_file, 'rb') as file:
    restored_params = pickle.load(file)


def pytree_compare(tree1, tree2):
    structure_equal = jax.tree_util.tree_flatten(tree1)[1] == jax.tree_util.tree_flatten(tree2)[1]
    entries_equal = np.all(jax.tree_util.tree_flatten(jax.tree_util.tree_map(np.array_equal, tree1, tree2))[0])
    return structure_equal and entries_equal


pytree_compare(params, restored_params)
"""


