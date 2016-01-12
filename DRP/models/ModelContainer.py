"""A module containing only the ModelContainer class."""
from django.db import models
from django.conf import settings
from django.db import transaction
from django.core.exceptions import ValidationError
from numpy import average
from itertools import chain
import random
import datetime
import importlib
import operator
import os
from rxnDescriptors import BoolRxnDescriptor, OrdRxnDescriptor, NumRxnDescriptor, CatRxnDescriptor
from StatsModel import StatsModel

visitorModules = {library:importlib.import_module(settings.STATS_MODEL_LIBS_DIR + "."+ library) for library in settings.STATS_MODEL_LIBS}

splitters = {splitter:importlib.import_module(settings.REACTION_DATASET_SPLITTERS_DIR + "." + splitter) for splitter in settings.REACTION_DATASET_SPLITTERS}
#TODO: set availability of manual splitting up

TOOL_CHOICES = tuple((key, tuple(tool for tool in library.tools)) for key,library in visitorModules.items())


class PredictsDescriptorsAttribute(object):

    def __get__(self, modelContainer, modelContainerType=None):
        return chain(modelContainer.predboolrxndescriptor_set.all(), modelContainer.predordrxndescriptor_set.all(), modelContainer.predcatrxndescriptor_set.all(), modelContainer.prednumrxndescriptor_set.all())

    def __set__(self, modelContainer, descriptors):
        modelContainer.predboolrxndescriptor_set.clear()
        modelContainer.predordrxndescriptor_set.clear()
        modelContainer.predcatrxndescriptor_set.clear()
        modelContainer.prednumrxndescriptor_set.clear()
        for descriptor in descriptors:
            descriptor.stats_modelContainer = modelContainer
            descriptors.save()

    def __delete__(self, modelContainer):
        modelContainer.predboolrxndescriptor_set.clear()
        modelContainer.predordrxndescriptor_set.clear()
        modelContainer.predcatrxndescriptor_set.clear()
        modelContainer.prednumrxndescriptor_set.clear()

class DescriptorAttribute(object):

    def __get__(self, modelContainer, modelContainerType=None):
        return chain(modelContainer.boolRxnDescriptors.all(), modelContainer.ordRxnDescriptors.all(), modelContainer.catRxnDescriptors.all(), modelContainer.numRxnDescriptors.all())

    def __set__(self, modelContainer, descriptors):
        modelContainer.boolRxnDescriptors.clear()
        modelContainer.ordRxnDescriptors.clear()
        modelContainer.catRxnDescriptors.clear()
        modelContainer.numRxnDescriptors.clear()
        for descriptor in descriptors:
            desc = None
            try:
                desc = BoolRxnDescriptor.objects.get(id=descriptor.id)
                modelContainer.boolRxnDescriptors.add(desc)
            except BoolRxnDescriptor.DoesNotExist:
                pass

            try:
                desc = OrdRxnDescriptor.objects.get(id=descriptor.id)
                modelContainer.ordRxnDescriptors.add(desc)
            except OrdRxnDescriptor.DoesNotExist:
                pass

            try:
                desc = CatRxnDescriptor.objects.get(id=descriptor.id)
                modelContainer.catRxnDescriptors.add(desc)
            except CatRxnDescriptor.DoesNotExist:
                pass

            try:
                desc = NumRxnDescriptor.objects.get(id=descriptor.id)
                modelContainer.numRxnDescriptors.add(desc)
            except NumRxnDescriptor.DoesNotExist:
                pass

            if desc is None:
                raise ValueError('An invalid object was assigned as a descriptor')

    def __delete__(self, modelContainer):
        modelContainer.boolRxnDescriptors.clear()
        modelContainer.numRxnDescriptors.clear()
        modelContainer.catRxnDescriptors.clear()
        modelContainer.ordRxnDescriptors.clear()

class OutcomeDescriptorAttribute(object):

    def __get__(self, modelContainer, modelContainerType=None):
        return chain(modelContainer.outcomeBoolRxnDescriptors.all(), modelContainer.outcomeOrdRxnDescriptors.all(), modelContainer.outcomeCatRxnDescriptors.all(), modelContainer.outcomeNumRxnDescriptors.all())

    def __set__(self, modelContainer, descriptors):
        modelContainer.outcomeBoolRxnDescriptors.clear()
        modelContainer.outcomeOrdRxnDescriptors.clear()
        modelContainer.outcomeCatRxnDescriptors.clear()
        modelContainer.outcomeNumRxnDescriptors.clear()
        for descriptor in descriptors:
            desc = None
            try:
                desc = BoolRxnDescriptor.objects.get(id=descriptor.id)
                modelContainer.outcomeBoolRxnDescriptors.add(desc)
            except BoolRxnDescriptor.DoesNotExist:
                pass

            try:
                desc = OrdRxnDescriptor.objects.get(id=descriptor.id)
                modelContainer.outcomeOrdRxnDescriptors.add(desc)
            except OrdRxnDescriptor.DoesNotExist:
                pass

            try:
                desc = CatRxnDescriptor.objects.get(id=descriptor.id)
                modelContainer.outcomeCatRxnDescriptors.add(desc)
            except CatRxnDescriptor.DoesNotExist:
                pass

            try:
                desc = NumRxnDescriptor.objects.get(id=descriptor.id)
                modelContainer.outcomeNumRxnDescriptors.add(desc)
            except NumRxnDescriptor.DoesNotExist:
                pass

            if desc is None:
                raise ValueError('An invalid object was assigned as a descriptor')

            pred_descriptor = desc.createPredictionDescriptor(modelContainer)
            pred_descriptor.save()

    def __delete__(self, modelContainer):
        modelContainer.outcomeBoolRxnDescriptors.clear()
        modelContainer.outcomeNumRxnDescriptors.clear()
        modelContainer.outcomeCatRxnDescriptors.clear()
        modelContainer.outcomeOrdRxnDescriptors.clear()

class ModelContainer(models.Model):

    """A class for describing a group of statistical models."""

    class Meta:
        app_label = 'DRP'


    description = models.TextField()
    active = models.BooleanField('Is this the active model?', default=False)
    library = models.CharField(
        max_length=200, choices=tuple((lib, lib) for lib in settings.STATS_MODEL_LIBS))
    tool = models.CharField(
        max_length=200, choices=tuple((tool, tool) for tool in TOOL_CHOICES))
    splitter = models.CharField(
        max_length=200, choices=tuple((splitter, splitter) for splitter in settings.REACTION_DATASET_SPLITTERS), blank=True, null=True)
    built = models.BooleanField('Has the build procedure been called with this container?', editable=False, default=False)

    descriptors = DescriptorAttribute()
    boolRxnDescriptors = models.ManyToManyField(BoolRxnDescriptor)
    ordRxnDescriptors = models.ManyToManyField(OrdRxnDescriptor)
    catRxnDescriptors = models.ManyToManyField(CatRxnDescriptor)
    numRxnDescriptors = models.ManyToManyField(NumRxnDescriptor)
    """The input descriptors for the model."""

    outcomeDescriptors = OutcomeDescriptorAttribute()
    outcomeBoolRxnDescriptors = models.ManyToManyField(BoolRxnDescriptor, related_name='outcomeForModels')
    outcomeOrdRxnDescriptors = models.ManyToManyField(OrdRxnDescriptor, related_name='outcomeForModels')
    outcomeCatRxnDescriptors = models.ManyToManyField(CatRxnDescriptor, related_name='outcomeForModels')
    outcomeNumRxnDescriptors = models.ManyToManyField(NumRxnDescriptor, related_name='outcomeForModels')
    """The descriptors which are being used as outcomes for this model.

    For models which make predictions about descriptors, it is probably
    most appropriate to make the descriptor label for the predicted
    descriptor related to the model, for instance for an outcome
    descriptor called "outcome", you might consider
    "outcome_predicted_by_model_id_1", where 1 is the model primary
    key.
    """

    predictsDescriptors = PredictsDescriptorsAttribute()
    """The descriptors which this model predicts values for."""


    fully_trained = models.ForeignKey("DRP.StatsModel", null=True)

    def __init__(self, library, tool, splitter=None, reactions=None, trainingSets=None, testSets=None):
        super(ModelContainer, self).__init__(splitter=splitter, library=library, tool=tool)


        self.reactions = reactions
        self.splitter = splitter
        self.reactions = reactions
        self.trainingSets = trainingSets
        self.testSets = testSets

    def clean(self):
        if self.tool not in visitorModules[self.library].tools:
            raise ValidationError('Selected tool does not exist in selected library', 'wrong_library')
        if getattr(visitorModules[self.library], self.tool).maxResponseCount is not None:
            if getattr(visitorModules[self.library], self.tool).maxResponseCount < self.outcomeDescriptors.count():
                raise ValidationError('Selected tool cannot accept this many responses, maximum is {}', 'too_many_responses', tuple(visitorModules[self.library], self.tool).maxResponseCount)
        if self.splitter is None ^ self.reactions is None:
            raise ValidationError('A full set of reactions must be supplied with a splitter', 'argument_mismatch')
        elif self.training is None:
            raise ValidationError('Either a splitter or a training set should be provided.', 'argument_mismatch')

    def build(self, predictors, response):

        if self.built:
            raise Exception("")

        self.descriptors = predictors
        self.outcomeDescriptors = response

        if self.splitter is not None:
            splitterObj = splitters[self.splitter].Splitter("{}_{}_{}".format(self.library, self.tool, self.pk))
            data_splits = splitterObj.split(self.reactions)

            self.trainingSets = [dataset_tuple[0] for dataset_tuple in data_splits]
            self.testSets = [dataset_tuple[1] for dataset_tuple in data_splits]

        resDict = {} #set up a prediction results dictionary. Hold on tight. This gets hairy real fast.

        for trainingSet, testSet in zip(self.trainingSets, self.testSets):
            statsModel = StatsModel(container=self, trainingSet=trainingSet)
            modelVisitor = getattr(visitorModules[self.library], self.tool)(statsModel)
            # Train the model.
            statsModel.startTime = datetime.datetime.now()
            fileName = os.path.join(settings.MODEL_DIR, '{}_{}'.format(self.pk, statsModel.pk))
            whitelist = [d.csvHeader for d in chain(self.descriptors, self.outcomeDescriptors)]
            modelVisitor.train(trainingSet.reactions.all(), whitelist, fileName)
            statsModel.fileName = fileName
            statsModel.endTime = datetime.datetime.now()
            statsModel.save()

            # Test the model.
            statsModel.testSets.add(testSet)
            predictions = modelVisitor.predict(testSet.reactions.all(), whitelist)
            newResDict = self._storePredictionComponents(predictions, statsModel)

            # Update the result-dictionary.
            for reaction, responseDict in newResDict:
                for response, outcome in responseDict:
                    if reaction not in resDict:
                        resDict[reaction] = {}
                    if response not in resDict[reaction]:
                        resDict[reaction][response] = {}

                    if outcome not in resDict[reaction][response]:
                        resDict[reaction][response][outcome] = 0
                    resDict[reaction][response][outcome] += 1


        self._storePredictions(resDict)
        self.built = True

    def _storePredictionComponents(self, predictions, statsModel):
        resDict = {}
        for response, outcome in predictions.items():
            for reaction, outcome in outcome:
                if reaction not in resDict:
                    resDict[reaction] = {}
                if response not in resDict[reaction]:
                    resDict[reaction][response] = {}
                if outcome not in resDict[reaction][response]:
                    resDict[reaction][response][outcome] = 0
                resDict[reaction][response][outcome] += 1

                predDesc = response.createPredictionDescriptor(self, statsModel)
                predDesc.save()

                val = predDesc.createValue(reaction, outcome)
                val.save()

        return resDict

    def _storePredictions(self, resDict):
        for reaction in resDict:
            for response in resDict[reaction]:
                predDesc = response.createPredictionDescriptor(self)
                predDesc.save()
                if isinstance(response, NumRxnDescriptor):
                    #estimate the weighted average of the estimates from the component models
                    values = tuple(value for value in resDict[reaction][response].keys())
                    weights = tuple(weight for value, weight in resDict[reaction][response].items())
                    val = predDesc.createValue(reaction, average(values, weights=weights))
                else:
                    #find the 'competitors' with the number of votes equal to the maximum and pick one at random for categorical variable types
                    votes = sorted(resDict[reaction][response].items(), key=operator.itemgetter(1), reverse=True)
                    maxVotes = max(v for k, v in votes)
                    winners = tuple(item for item in votes if item[1] == maxVotes)
                    winner = random.choice(winners)
                    val = predDesc.createValue(reaction, winner)
                val.save()

    def predict(self, reactions):
        if self.built:
            for model in self.statsmodel_set:
                modelVisitor = getattr(visitorModules[self.library], self.tool)(model)
                predictions = modelVisitor.predict(reactions)
                resDict = self._storePredictionComponents(predictions, model)
                self._storePredictions(resDict)
        else:
            raise RuntimeError('A model container cannot be used to make predictions before the build method has been called')

    @transaction.atomic
    def save(self, *args, **kwargs):
        super(ModelContainer, self).save(*args, **kwargs)

    def summarize(self):
        """CAUTION: This is a temporary development function. Do not rely on it. """
        """Return a string containing the Confusion Matrices for all stats_models."""
        summaries = "\nK-Fold Validation:\n"
        for model in self.statsmodel_set.all():
            for d in model.predictsDescriptors:
                summaries += "{}\n".format(d.summarize(model))
        return summaries
