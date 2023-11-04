import numpy as np

from sklearn import datasets
from sklearn.model_selection import train_test_split

RND_SEED = None

def load_data():
    iris = datasets.load_iris()
    X = iris.data[:,[2,3]] # only two features
    y = iris.target
    return X, y


X, y = load_data()
X_train, X_test, y_train, y_test = train_test_split(
    X, y, test_size=0.3, random_state=RND_SEED, stratify=y
)
# stratify: training and test subsets have the same proportions of class
# labels as the input dataset. 

print("Labels counts in y: ", np.bincount(y))
print("Labels counts in y_train = ", np.bincount(y_train))
print("Labels counts in y_test = ", np.bincount(y_test))

from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.pipeline import make_pipeline


model = make_pipeline(
    StandardScaler(),
    LogisticRegression(solver="lbfgs", multi_class="ovr", random_state=RND_SEED)
)
model.fit(X_train, y_train)

y_pred = model.predict(X_test)
y_pred2 = model.predict(X_train)
print("Misclassified examples in test data: %d" % (y_test != y_pred).sum())
print("Misclassified examples in train data: %d" % (y_train != y_pred2).sum())

#from sklearn.metrics import accuracy_score
#print("Accuracy = %.3f" % accuracy_score(y_test, y_pred))
# Alternative
print("Accuracy (train data) = %.3f" % model.score(X_train, y_train))
print("Accuracy (test data) = %.3f" % model.score(X_test, y_test))

from matplotlib.colors import ListedColormap
import matplotlib.pyplot as plt
def plot_decision_regions(X, y, classifier, test_idx=None,
                          resolution=0.02):
    # setup marker generator and color map
    markers = ("o", "s", "^", "v", "<")
    colors = ("red", "blue", "lightgreen", "gray", "cyan")
    cmap = ListedColormap(colors[:len(np.unique(y))])
    # plot the decision surface
    x1_min, x1_max = X[:, 0].min() - 1, X[:, 0].max() + 1
    x2_min, x2_max = X[:, 1].min() - 1, X[:, 1].max() + 1
    xx1, xx2 = np.meshgrid(np.arange(x1_min, x1_max, resolution),
                           np.arange(x2_min, x2_max, resolution))
    lab = classifier.predict(np.array([xx1.ravel(), xx2.ravel()]).T)
    lab = lab.reshape(xx1.shape)
    plt.contourf(xx1, xx2, lab, alpha=0.3, cmap=cmap)
    plt.xlim(xx1.min(), xx1.max())
    plt.ylim(xx2.min(), xx2.max())
    # plot class examples
    for idx, cl in enumerate(np.unique(y)):
        plt.scatter(x=X[y == cl, 0],
                    y=X[y == cl, 1],
                    alpha=0.8,
                    c=colors[idx],
                    marker=markers[idx],
                    label=f"Class {cl}",
                    edgecolor="black")
    # highlight test examples
    if test_idx:
        # plot all examples
        X_test, y_test = X[test_idx, :], y[test_idx]
        plt.scatter(X_test[:, 0], X_test[:, 1],
                    c="none", edgecolor="black", alpha=1.0,
                    linewidth=1, marker="o",
                    s=100, label="Test set")



X_combined = np.vstack((X_train, X_test))
y_combined = np.hstack((y_train, y_test))
plot_decision_regions(
    X=X_train, y=y_train,
    classifier=model
)
plt.xlabel("Petal length")
plt.ylabel("Petal width")
plt.legend(loc="upper left")
plt.tight_layout()
plt.show()