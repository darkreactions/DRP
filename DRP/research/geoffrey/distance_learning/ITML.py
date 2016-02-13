#!/usr/bin/env python
import numpy as np
from metric_learn import ITML
from sklearn.datasets import load_iris
from DRP.models import PerformedReaction, ModelContainer, Descriptor, rxnDescriptorValues
from DRP.models.rxnDescriptorValues import OrdRxnDescriptorValue, NumRxnDescriptorValue, BoolRxnDescriptorValue, CatRxnDescriptorValue
from django.db.models import Q
import operator
from sys import argv

class ITML(AbstractDistanceLearner):
  def __init__(self, *args, **kwargs):
    super(self.__class__, self).__init__(*args, **kwargs)
    
    def _train(predictors, responses):
      iris_data = load_iris()
      X = iris_data['data']
      Y = iris_data['target']

      itml = ITML()

      num_constraints = 200
      C = ITML.prepare_constraints(Y, X.shape[0], num_constraints)
      itml.fit(X, C, verbose=False)
