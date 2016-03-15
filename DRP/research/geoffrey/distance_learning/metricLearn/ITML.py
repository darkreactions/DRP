#!/usr/bin/env python
from metric_learn import ITML as ml_ITML
import numpy as np
from sklearn.metrics import pairwise_distances
from DRP.research.geoffrey.distance_learning.metricLearn.AbstractMetricLearnDistanceLearner import AbstractMetricLearnDistanceLearner, logger


class MetricVisitor(AbstractMetricLearnDistanceLearner):
    def __init__(self, num_constraints, *args, **kwargs):
        self.metric_object = ml_ITML()
        self.num_constraints = num_constraints if num_constraints is not None else 200
        super(MetricVisitor, self).__init__(*args, **kwargs)


    
    def train(self, reactions, predictor_headers, response_headers, filename):
        print "Preparing arrays"
        data, labels = self._prepareArrays(reactions, predictor_headers, response_headers)
        old_settings = np.seterr(divide='raise') # we don't want division by zero to pass

        # This is how metric learn determines bounds internally
        # but the lower bound can be zero this way (especially for low-dimensional data)
        # which causes divide by zero errors
        print "Calculating bounds"
        pair_dists = pairwise_distances(data)
        bounds = np.percentile(pair_dists, (5, 95))
        # the extra check ensures against divide-by-zero errors later
        if bounds[0] == 0:
            bounds[0] = min(pair_dists[np.nonzero(pair_dists)])
            print "Lowerbound was 0. Set to {}".format(bounds[0])
        
        print "Preparing {} constraints with bounds of ({}, {})".format(self.num_constraints, bounds[0], bounds[1])
        constraints = self.metric_object.prepare_constraints(labels, data.shape[0], self.num_constraints)
        print "Fitting"
        self.metric_object.fit(data, constraints, bounds=bounds)
        
        self.save(filename)
        np.seterr(**old_settings)
        
        print "Transforming training set"
        return self.metric_object.transform()


            
