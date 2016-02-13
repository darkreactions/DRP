#!/usr/bin/env python
import numpy as np
from metric_learn import ITML as ml_ITML
from DRP.models import PerformedReaction, Descriptor, rxnDescriptorValues
from DRP.research.geoffrey.distance_learning.AbstractDistanceLearner import AbstractDistanceLearner, logger

class ITML(AbstractDistanceLearner):
    def __init__(self, *args, **kwargs):
        super(self.__class__, self).__init__(*args, **kwargs)
    
    def train(self, reactions, predictorHeaders, responseHeaders):
        X = reactions.toNPArray(expanded=True, whitelistHeaders=predictorHeaders)
        exit(1)
        #iris_data = load_iris()
        #X = iris_data['data']
        #Y = iris_data['target']
        
        #itml = ITML()
        
        #num_constraints = 200
        #C = ITML.prepare_constraints(Y, X.shape[0], num_constraints)
        #itml.fit(X, C, verbose=False)
    
