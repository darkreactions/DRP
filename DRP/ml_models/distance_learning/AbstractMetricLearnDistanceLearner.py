"""Placeholder summary line.

Placeholder description.
"""

from DRP.research.geoffrey.distance_learning.AbstractDistanceLearner import AbstractDistanceLearner, logger
import numpy as np
from sklearn.preprocessing import Imputer
from cPickle import dump, load


class AbstractMetricLearnDistanceLearner(AbstractDistanceLearner):
    """Place holder doc string."""

    maxResponseCount = 1

    def _prepareArrays(self, reactions, predictor_headers, response_headers):
        data = reactions.toNPArray(
            expanded=True, whitelistHeaders=predictor_headers, missing=np.nan)
        labels = reactions.toNPArray(
            expanded=True, whitelistHeaders=response_headers).flatten()

        data = Imputer(copy=False).fit_transform(data)

        try:
            assert(data.dtype == np.float64)
        except AssertionError:
            raise TypeError("Data is not of the type float.")

        return data, labels

    def transform(self, reactions, predictor_headers):
        """Place holder doc string."""
        data = self._prepareArrays(reactions, predictor_headers, [])[0]
        transformed = self.metric_object.transform(data)
        return transformed

    def metric(self):
        """Place holder doc string."""
        return self.metric_object.metric()

    def dist(self, x, y):
        """Place holder doc string."""
        x = np.array(x)
        y = np.array(y)
        Mah = self.metric_object.metric()
        dif = x - y

        return dif.T.dot(Mah).dot(dif)

    def save(self, fileName):
        """Place holder doc string."""
        with open(fileName, 'wb') as f:
            dump(self.metric_object, f)

    def recover(self, fileName):
        """Place holder doc string."""
        with open(fileName, 'rb') as f:
            self.metric_object = load(f)
