from abc import ABCMeta, abstractmethod
import os
from DRP.models import StatsModel, PerformedReaction, TrainingSet, TestSet, TestSetRelation, Descriptor
from DRP.models.predRxnDescriptors import PredBoolRxnDescriptor, PredOrdRxnDescriptor, PredNumRxnDescriptor, PredCatRxnDescriptor
from DRP.models.rxnDescriptors import BoolRxnDescriptor, OrdRxnDescriptor, NumRxnDescriptor, CatRxnDescriptor
from DRP.models.rxnDescriptorValues import BoolRxnDescriptorValue, OrdRxnDescriptorValue, NumRxnDescriptorValue, CatRxnDescriptorValue
from django.conf import settings
from django.core.files import File
import logging
#TODO: set attribute methods to be transactions
#TODO: set descriptors to forbid the word predicted
#TODO: input logging options into the settings files

logger = logging.getLogger(__name__)

class PredictorsAttribute(object):

    def __get__(self, visitor, visitorType=None):
        return visitor.stats_model.descriptors

    def __set__(self, visitor, descriptors)
        descriptors = list(descriptors)
        if len(descriptors) < 1:
            raise ValueError('At least one descriptor must be set')
        visitor.stats_model.descriptors = descriptors
        visitor.stats_model.save()

    def __delete__(self, visitor):
        del visitor.stats_model.descriptors
        visitor.stats_model.save()

class ResponsesAttribute(object):

    def __get__(self, visitor, visitorType=None):
        return self.stats_model.outcomeDescriptors 

    def __set__(self, visitor, descriptors)
        descriptors = list(descriptors)
        if len(descriptors) < 1:
            raise ValueError('Response descriptors cannot be empty!')
        visitor.stats_model.outcomeDescriptors = descriptors
        visitor.stats_model.save()

    def __delete__(self, visitor):
        del visitor.stats_model.descriptors
        visitor.stats_model.save()

class TrainingDataAttribute(object):
    
    def __get__(self, visitor, visitorType=None):
        return visitor.stats_model.trainingSet

    def __set__(self, visitor, dataSet):
        visitor.stats_model.trainingSet = dataSet
        visitor.stats_model.save()

class TestingDataAttribute(object):

    def __get__(self, visitor, visitorType=None):
        return visitor.stats_model.testSets

    def __set__(self, visitor, dataSets):
        visitor.stats_model.testSets.clear()
        for dataSet in dataSets:
            visitor.stats_model.testSets.add(dataset)
        visitor.stats_model.save()

    def __delete__(self, visitor)
        visitor.stats_model.testSets.clear()
        visitor.stats_model.save()


class AbstractModelVisitor(object):
    __metaclass__ = ABCMeta

    maxResponseCount = None

    def __init__(self, statsModel):
        self.statsModel = statsModel

    @abstractmethod
    def train(self, reactions, descriptorHeaders, filePath):
        """A function meant to be overridden by actual ModelVisitor classes.
              The `_train` method should prepare the machine learning model for
              classification and save that model if necessary."""

    @abstractmethod
    def predict(self, reactions, descriptorHeaders):
        """Return a dictionary where the key is the response descriptor being predicted
              and the value is an ordered list of predictions for that response where the
              ith prediction corresponds to the ith reaction.

              EG: {<NumRxnDescriptor> "outcome" }:[1,2,1,1]}"""
