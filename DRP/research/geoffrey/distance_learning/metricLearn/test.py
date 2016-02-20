from modshogun import LMNN as shogun_LMNN
from modshogun import RealFeatures, MulticlassLabels
import numpy as np
from metric_learn import LMNN
from sklearn.datasets import load_iris

iris_data = load_iris()
X = iris_data['data']
Y = iris_data['target']

lmnn = LMNN(k=5, learn_rate=1e-6)
lmnn.fit(X, Y, verbose=False)
