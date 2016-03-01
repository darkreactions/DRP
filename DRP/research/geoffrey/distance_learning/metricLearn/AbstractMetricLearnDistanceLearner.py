from DRP.research.geoffrey.distance_learning.AbstractDistanceLearner import AbstractDistanceLearner, logger
import numpy as np
from sklearn.preprocessing import Imputer
from cPickle import dump, load

class AbstractMetricLearnDistanceLearner(AbstractDistanceLearner):
    maxResponseCount = 1

    def _prepareArrays(self, reactions, predictor_headers, response_headers):
        data = reactions.toNPArray(expanded=True, whitelistHeaders=predictor_headers, missing=np.nan)
        labels = reactions.toNPArray(expanded=True, whitelistHeaders=response_headers).flatten()

        data = Imputer(copy=False).fit_transform(data)

        try:
            assert(data.dtype == np.float64)
        except AssertionError:
            raise TypeError("Data is not of the type float.")

        return data, labels

    def transform(self, reactions, predictor_headers):
        data = self._prepareArrays(reactions, predictor_headers, [])[0]
        transformed = self.metric_object.transform(data)
        return transformed

    def metric(self):
        return self.metric_object.metric()

    def dist(self, x, y):
        x = np.array(x)
        y = np.array(y)
        Mah = self.metric_object.metric()
        dif = x - y

        return dif.T.dot(Mah).dot(dif)
        
    def save(self, fileName):
        with open(fileName, 'wb') as f:
            dump(self.metric_object, f)

    def recover(self, fileName):
        with open(fileName, 'rb') as f:
            self.metric_object = load(f)
