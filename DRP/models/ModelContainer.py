"""A module containing only the ModelContainer class."""
from django.db import models
from django.conf import settings
from django.db import transaction
from django.core.exceptions import ValidationError
import datetime
import importlib
visitorModules = {library:importlib.import_module(library) for library in settings.STATS_MODEL_LIBS}
splitters = {splitterName:importlib.import_module(splitter) for splitter in settings.REACTION_DATASET_SPLITTERS}
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
            if isinstance(descriptor, BoolRxnDescriptor):
                modelContainer.boolRxnDescriptors.add(descriptor):
            elif isinstance(descriptor, OrdRxnDescriptor):
                modelContainer.ordRxnDescriptors.add(descriptor)
            elif isinstance(descriptor, CatRxnDescriptor):
                modelContainer.catRxnDescriptors.add(descriptor)
            elif isinstance(descriptor, NumRxnDescriptor):
                modelContainer.numRxnDescriptors.add(descriptor)
            else:
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
            if isinstance(descriptor, BoolRxnDescriptor):
                modelContainer.outcomeBoolRxnDescriptors.add(descriptor):
                pred_descriptor = PredBoolRxnDescriptor()
            elif isinstance(descriptor, OrdRxnDescriptor):
                modelContainer.outcomeOrdRxnDescriptors.add(descriptor)
                pred_descriptor = PredOrdRxnDescriptor()
                pred_descriptor.maximum = descriptor.maximum
                pred_descriptor.minimum = descriptor.minimum
            elif isinstance(descriptor, CatRxnDescriptor):
                pred_descriptor = PredCatRxnDescriptor()
                modelContainer.outcomeCatRxnDescriptors.add(CatRxnDescriptor)
            elif isinstance(descriptor, NumRxnDescriptor):
                modelContainer.outcomeNumRxnDescriptors.add(NumRxnDescriptor)
                pred_descriptor = PredNumRxnDescriptor()
                pred_descriptor.maximum = descriptor.maximum
                pred_descriptor.minimum = descriptor.minimum
            else:
                raise ValueError('An invalid object was assigned as a descriptor')
            pred_descriptor.heading = descriptor.heading + modelContainer.predictionSuffix
            pred_descriptor.name = descriptor.name + modelContainer.predictionSuffix

            pred_descriptor.prediction_of = descriptor
            pred_descriptor.stats_modelContainer = modelContainer
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
        max_length=200, choices=settings.STATS_MODEL_LIBS)
    tool = models.CharField(
        max_length=200, choices=TOOL_CHOICES)
    splitter = models.CharField(
        max_length=200, choices=settings.REACTION_DATASET_SPLITTERS, blank=True, null=True)
    built = models.Booleanfield('Has the build procedure been called with this container?', editable=False, default=False)

    descriptors = DescriptorAttribute()
    boolRxnDescriptors = models.ManyToManyField(BoolRxnDescriptor)
    ordRxnDescriptors = models.ManyToManyField(OrdRxnDescriptor)
    catRxnDescriptors = models.ManyToManyField(CatRxnDescriptor)
    numRxnDescriptors = models.ManyToManyField(NumRxnDescriptor)
    """The input descriptors for the model."""

    outcomeDescriptors = OutComeDescriptorAttribute()
    outcomeBoolRxnDescriptors = models.ManyToManyField(BoolRxnDescriptor, related_name='outcomeForModels'))
    outcomeOrdRxnDescriptors = models.ManyToManyField(OrdRxnDescriptor, related_name='outcomeForModels'))
    outcomeCatRxnDescriptors = models.ManyToManyField(CatRxnDescriptor, related_name='outcomeForModels'))
    outcomeNumRxnDescriptors = models.ManyToManyField(NumRxnDescriptor, related_name='outcomeForModels'))
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

    def __init__(self, responses, splitter, library, tool, splitter=None, reactions=None, trainingSets=None, testSets=None):
        super(ModelContainer, self).__init__(splitter=splitter, library=library, tool=tool)
        self.reactions = reactions
        self.responses = responses
        self.splitter = splitter
        self.reactions = reactions
        self.trainingSets = trainingSets
        self.testSets = testSets

    def clean(self):
        if self.tool not in visitorModules[self.library].tools:
            raise ValidationError('Selected tool does not exist in selected library', 'wrong_library')
        if self.splitter is None ^ self.reactions is None:
            raise ValidationError('A full set of reactions must be supplied with a splitter', 'argument_mismatch')
        elif self.training is None
            raise ValidationError('Either a splitter or a training set should be provided.', 'argument_mismatch') 

    def build(self):
        if not self.built:
            if self.splitter is not None:
                data = splitters[self.splitter].Splitter("{}_{}_{}".format(self.library, self.tool, self.pk)).split(self.reactions)
                self.trainingSets = list(d[0] for d in data)
                self.testSets = list(d[1] for d in data)
            for trainingSet, testSets in zip(self.trainingSets, selt.testSets):
                statsModel = StatsModel(container=self, trainingSet=trainingSet)
                modelVisitor = visitorModules[self.library].ModelVisitor(statsModel)
                statsModel.startTime = datetime.datetime.now()
                fileName = os.path.join(settings.MODEL_DIR, '{}_{}'.format(self.pk, stats_model.pk))
                whitelist = [d.csvHeader for d in chain(self.descriptors, self.outcomeDescriptors)]
                modelVisitor.train(trainingSet.reactions.all(), whiteList, fileName)
                self.stats_model.fileName = fileName
                statsModel.endTime = datetime.datetime.now()
                statsModel.save()
                for testSet in testSets:
                    statsModel.testSets.add(testSet)  
                    modelVisitor.predict(testSet.reactions.all()) #TODO: START HERE
            self.built = True

    def test(self, testSet):
        

    @transaction.atomic
    def save(self, *args, **kwargs):
        super(ModelContainer, self).save(*args, **kwargs)
            modelVisitor.test()

    def summarize(self):
        """CAUTION: This is a temporary development function. Do not rely on it. """
        """Return a string containing the Confusion Matrices for all stats_models."""
        summaries = "\nK-Fold Validation:\n"
        for model in self.statsmodel_set.all():
            for d in model.predictsDescriptors:
                summaries += "{}\n".format(d.summarize(model))

        return summaries

    def duplicate(self):
