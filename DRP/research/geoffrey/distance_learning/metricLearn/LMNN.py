#!/usr/bin/env python
from metric_learn import LMNN as ml_LMNN
from DRP.research.geoffrey.distance_learning.metricLearn.AbstractMetricLearnDistanceLearner import AbstractMetricLearnDistanceLearner, logger

class LMNN(AbstractMetricLearnDistanceLearner):
    def __init__(self, *args, **kwargs):
        self.metric_object = ml_LMNN()
        super(LMNN, self).__init__(*args, **kwargs)
    
    def train(self):
        self._prepareArrays()
        self.metric_object.fit(self.data, self.labels)
