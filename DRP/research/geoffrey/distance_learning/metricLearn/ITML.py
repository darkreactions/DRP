#!/usr/bin/env python
import numpy as np
from metric_learn import ITML as ml_ITML
from sklearn.preprocessing import Imputer
from sklearn.metrics import pairwise_distances
from DRP.models import PerformedReaction, Descriptor, rxnDescriptorValues
from DRP.research.geoffrey.distance_learning.AbstractDistanceLearner import AbstractDistanceLearner, logger
import cPickle as pickle

class ITML(AbstractDistanceLearner):
    def __init__(self, *args, **kwargs):
        super(self.__class__, self).__init__(*args, **kwargs)
        self.matrix = None
    
    def train(self, reactions, predictorHeaders, responseHeaders):
        data = reactions.toNPArray(expanded=True, whitelistHeaders=predictorHeaders, missing=np.nan)
        # TODO XXX: unset missing and ensure that all labels are defined.
        # We shouldn't be using reactions with undefined labels anyway
        labels = reactions.toNPArray(expanded=True, whitelistHeaders=responseHeaders, missing=False).flatten()

        data = Imputer().fit_transform(data)

        bounds = np.percentile(pairwise_distances(data),(5,95))
        # this is how metric-learn sets bounds internally if none are given.
        # determine them here to make sure they're not 0!

        bounds[0] = max(bounds[0], 0.1)

        np.seterr(all='raise') # we don't want division by zero to pass
        try:
            assert(data.dtype == np.float64)
        except AssertionError:
            raise TypeError("Data is not of the type float.")

        itml = ml_ITML()

        num_constraints = 10000
        constraints = ml_ITML.prepare_constraints(labels, data.shape[0], num_constraints)
        itml.fit(data, constraints, bounds=bounds)

        self.matrix = itml.metric()
        
    def save(self, writeable):
        pickle.dump(self, writeable)
    
