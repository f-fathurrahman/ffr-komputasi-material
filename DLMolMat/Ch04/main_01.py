import pandas as pd
import mordred, mordred.descriptors
import rdkit

import numpy as np
import jax
import jax.numpy as jnp

toxdata = pd.read_csv("../DATASET/clintox.csv")
print(toxdata.head())

calc = mordred.Calculator(mordred.descriptors, ignore_3D=True)
# Samples
molecules = [rdkit.Chem.MolFromSmiles(smi) for smi in toxdata.smiles]

# the invalid molecules were None, so we'll just
# use the fact the None is False in Python
valid_mol_idx = [bool(m) for m in molecules]
valid_mols = [m for m in molecules if m]
print("Number of valid molecules: ", len(valid_mols))

# Requires pandas 1.x
#print("Calculating features:") # This takes quite long time
#features = calc.pandas(valid_mols, nproc=1)
#print("Finished calculating features")

#features = pd.read_csv("../DATASET/features_mordred.csv", low_memory=True)
features = pd.read_pickle("../DATASET/features_mordred.pkl")

features.dropna(inplace=True, axis=1)

labels = toxdata[valid_mol_idx].FDA_APPROVED
features -= features.mean(numeric_only=True)
features /= features.std(numeric_only=True)

# we have some nans in features, likely because std was 0
features.dropna(inplace=True, axis=1)
print(f"We have {len(features.columns)} features per molecule")


def perceptron(x, w, b):
    v = jnp.dot(x, w) + b
    y = jnp.where(v > 0, x=jnp.ones_like(v), y=jnp.zeros_like(v))
    return y



def loss(y, yhat):
    return jnp.mean(jnp.abs(y - yhat))


def loss_wrapper(w, b, x, y):
    yhat = perceptron(x, w, b)
    return loss(y, yhat)


loss_grad = jax.grad(loss_wrapper, (0, 1))


batch_size = 32
train_N = int(len(labels) * 0.8)


N = len(labels)
batch_idx = range(0, train_N, batch_size)
w = np.random.normal(size=len(features.columns))
b = 0.0

loss_grad = jax.grad(loss_wrapper, (0, 1))


test_x = features[train_N:].values.astype(np.float32)
test_y = labels[train_N:].values


loss_grad(w, b, test_x, test_y)


def bin_classifier(x, w, b):
    v = jnp.dot(x, w) + b
    y = jax.nn.sigmoid(v)
    return y


def cross_ent(y, yhat):
    return jnp.mean(-(y * jnp.log(yhat + 1e-10) + (1 - y) * jnp.log(1 - yhat + 1e-10)))


def loss_wrapper(w, b, x, y):
    yhat = bin_classifier(x, w, b)
    return cross_ent(y, yhat)


loss_grad = jax.grad(loss_wrapper, (0, 1))
w = np.random.normal(scale=0.01, size=len(features.columns))
b = 1.0

loss_progress = []
test_loss_progress = []
eta = 0.2
for epoch in range(5):
    for i in range(len(batch_idx) - 1):
        x = features[batch_idx[i] : batch_idx[i + 1]].values.astype(np.float32)
        y = labels[batch_idx[i] : batch_idx[i + 1]].values
        grad = loss_grad(w, b, x, y)
        w -= eta * grad[0]
        b -= eta * grad[1]
        loss_progress.append(loss_wrapper(w, b, x, y))
        test_loss_progress.append(loss_wrapper(w, b, test_x, test_y))

import matplotlib.pyplot as plt

plt.clf()
plt.plot(loss_progress, label="Training Loss")
plt.plot(test_loss_progress, label="Testing Loss")
plt.xlabel("Step")
plt.legend()
plt.ylabel("Loss")
plt.show()
plt.savefig("IMG_training_01.png")


# Metrics

def accuracy(y, yhat):
    # convert from prob to hard class
    hard_yhat = np.where(yhat > 0.5, np.ones_like(yhat), np.zeros_like(yhat))
    disagree = np.sum(np.abs(y - yhat))
    return 1 - disagree / len(y)

score_acc = accuracy(test_y, bin_classifier(test_x, w, b))
print("accuracy = ", score_acc)


# Alternative classifier
def alt_classifier(x):
    return np.ones((x.shape[0]))


score_acc_alt = accuracy(test_y, alt_classifier(test_x))
print("accuracy alt classifier = ", score_acc_alt)

def error_types(y, yhat, threshold):
    hard_yhat = np.where(yhat >= threshold, np.ones_like(yhat), np.zeros_like(yhat))
    # predicted 1, actually was 0 -> 1 (bool to remove predicted 0, actually was 1)
    fp = np.sum((hard_yhat - y) > 0)
    # predicted 0, actually was 1 -> 1 (bool to remove predicted 1, actually was 0)
    fn = np.sum((y - hard_yhat) > 0)
    return fp, fn

print("Alt Classifier (FP, FN) = ", error_types(test_y, alt_classifier(test_x), 0.5))
print("Trained Classifier (FP, FN) = ", error_types(test_y, bin_classifier(test_x, w, b), 0.5))

# Now we have a better sense of how our model does in comparison.
# The number of errors is indeed larger for our trained model, but
# it has a bit of balance between the two errors. What is more important?
# In our case, I would argue doing clinical trials that fail is worse than mistakenly
# not starting them. That is, false positives are worse than false negatives.
# Let's see if we can tune our threshold value to minimize false positives.

# Let's see if we can tune our threshold value to minimize false positives.
print("Threshold 0.7", error_types(test_y, bin_classifier(test_x, w, b), 0.7))
print("Threshold 0.9", error_types(test_y, bin_classifier(test_x, w, b), 0.9))
print("Threshold 0.95", error_types(test_y, bin_classifier(test_x, w, b), 0.95))
print("Threshold 0.99", error_types(test_y, bin_classifier(test_x, w, b), 0.99))


# By adjusting the threshold, we can achieve a balance of error more like what we desire for
# our model. We’re able to have 1 false positives in fact, at the cost of missing 218 of
# the molecules. Now are we still predicting positives? Are we actually going to get some
# true positives? We can measure that as well

total_pos = np.sum(test_y)
print(
    "Total positives:",
    total_pos,
    "Predicted Positives:",
    np.sum(bin_classifier(test_x, w, b) > 0.99),
)

# Yes, our model is actually capable of predicting if molecules will pass FDA clinical
# trials with as few false positives as possible (1). A model that is capable of this
# tuning is an example of a good model. Our other model, that predicts 1s, has good
# accuracy but we cannot adjust it or try to better balance type I and type II errors.


# ROC

unique_threshes = np.unique(bin_classifier(test_x, w, b))
fp = []
tp = []
total_pos = np.sum(test_y)
for ut in list(unique_threshes) + [-0.1, 1.01]:
    errors = error_types(test_y, bin_classifier(test_x, w, b), ut)
    fp.append(errors[0])
    tp.append(total_pos - errors[1])

# sort them so can plot as a line
idx = np.argsort(fp)
fpr = np.array(fp)[idx] / (len(test_y) - np.sum(test_y))
tpr = np.array(tp)[idx] / np.sum(test_y)

# now remove duplicate x-values
fpr_nd = []
tpr_nd = []
last = None
for f, t in zip(fpr, tpr):
    if last is None or f != last:
        last = f
        fpr_nd.append(f)
        tpr_nd.append(t)

plt.plot(fpr_nd, tpr_nd, "-o", label="Trained Model")
plt.plot([0, 1], [0, 1], label="Naive Classifier")
plt.ylabel("True Positive Rate")
plt.xlabel("False Positive Rate")
plt.legend()
plt.show()
plt.savefig("IMG_ROC_01.png")


# This plot nicely shows how our trained model is actually sensitive to threshold, so that we
# could choose to more carefully screen for false negative or false positives. The best curves
# fall to the top-left of this plot. Our naive classifier is where we return a fixed percentage
# of examples randomly as positive or negative. You can plot the area under this curve with an
# integration and this is a good way to measure classifier performance and correctly capture the
# effect of both false negatives and false positives. The area under the ROC curve is known as
# the ROC AUC score and is preferred to accuracy because it captures the balance of Type I and II errors.


