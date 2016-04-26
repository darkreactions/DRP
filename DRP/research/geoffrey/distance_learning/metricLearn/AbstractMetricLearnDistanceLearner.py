from DRP.research.geoffrey.distance_learning.AbstractDistanceLearner import AbstractDistanceLearner, logger
import numpy as np
from sklearn.preprocessing import Imputer


class AbstractMetricLearnDistanceLearner(AbstractDistanceLearner):
    maxResponseCount = 1
    metric_object = None

    def _prepareArrays(self):
        self.data = self.reactions.toNPArray(expanded=True, whitelistHeaders=self.predictorHeaders, missing=np.nan)
        # TODO XXX: unset missing and ensure that all labels are defined.
        # We shouldn't be using reactions with undefined labels anyway
        self.labels = self.reactions.toNPArray(expanded=True, whitelistHeaders=self.responseHeaders, missing=False).flatten()

        Imputer(copy=False).fit_transform(self.data)

        try:
            assert(self.data.dtype == np.float64)
        except AssertionError:
            raise TypeError("Data is not of the type float.")

    def transformer(self):
        return self.metric_object.transformer()

    def transform(self, X=None):
        return self.metric_object.transform(X)

    def metric(self):
        return self.metric_object.metric()

    def dist(self, x, y):
        x = np.array(x)
        y = np.array(y)
        Mah = self.metric_object.metric()
        dif = x - y

        return dif.T.dot(Mah).dot(dif)
        

