#!/usr/bin/env python
import numpy as np
from metric_learn import ITML as ml_ITML
from sklearn.preprocessing import Imputer
from sklearn.metrics import pairwise_distances
from DRP.models import PerformedReaction, Descriptor, rxnDescriptorValues
from DRP.research.geoffrey.distance_learning.AbstractDistanceLearner import AbstractDistanceLearner, logger

class ITML(AbstractDistanceLearner):
    def __init__(self, *args, **kwargs):
        super(self.__class__, self).__init__(*args, **kwargs)
    
    def train(self, reactions, predictorHeaders, responseHeaders):
        data = reactions.toNPArray(expanded=True, whitelistHeaders=predictorHeaders, missing=np.nan)
        labels = reactions.toNPArray(expanded=True, whitelistHeaders=responseHeaders, missing=False).flatten()

        imp = Imputer()
        data = imp.fit_transform(data)

        bounds = np.percentile(pairwise_distances(data),(5,95))
        # this is how metric-learn sets bounds internally if none are given.
        # determine them here to make sure they're not 0!

        bounds[0] = max(bounds[0], 0.1)

        np.seterr(all='raise') # we don't want division by zero to pass
        try:
            assert(data.dtype == np.float64)
        except AssertionError:
            raise TypeError("Data is not of the type float.")
        try:
            assert(labels.dtype == np.bool)
        except AssertionError:
            raise TypeError("Labels are not of the type bool.")

        itml = ml_ITML()

        num_constraints = 200
        print "generating constraints...",
        constraints = ml_ITML.prepare_constraints(labels, data.shape[0], num_constraints)
        print "constraints generated"
        itml.fit(data, constraints, bounds=bounds)
    
