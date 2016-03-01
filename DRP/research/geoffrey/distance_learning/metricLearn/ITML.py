#!/usr/bin/env python
from metric_learn import ITML as ml_ITML
import numpy as np
from sklearn.metrics import pairwise_distances
from DRP.research.geoffrey.distance_learning.metricLearn.AbstractMetricLearnDistanceLearner import AbstractMetricLearnDistanceLearner, logger
from cPickle import dump, load

class ITML(AbstractMetricLearnDistanceLearner):
    def __init__(self, *args, **kwargs):
        self.metric_object = ml_ITML()
        super(self.__class__, self).__init__(*args, **kwargs)


    
    def train(self, reactions, predictor_headers, response_headers, filename):
        data, labels = self._prepareArrays(reactions, predictor_headers, response_headers)
        num_constraints = num_constraints #record so we know what it was trained with
        old_settings = np.seterr(divide='raise') # we don't want division by zero to pass

        # This is how metric learn determines bounds internally
        # but the lower bound can be zero this way (especially for low-dimensional data)
        # which causes divide by zero errors
        pair_dists = pairwise_distances(data)
        bounds = np.percentile(pair_dists, (5, 95))
        # the extra check ensures against divide-by-zero errors later
        if bounds[0] == 0:
            bounds[0] = min(pair_dists[np.nonzero(pair_dists)])

        constraints = self.metric_object.prepare_constraints(labels, data.shape[0], num_constraints)
        self.metric_object.fit(data, constraints, bounds=bounds)

        self.save(filename)
        np.seterr(**old_settings)


            
