"""A module containing only the ModelContainer class."""
from django.db import models
from django.conf import settings
from django.db import transaction
from django.core.exceptions import ValidationError
import importlib
visitorModules = {library:importlib.import_module(library) for library in settings.STATS_MODEL_LIBS}
splitters = {splitterName:importlib.import_module(splitter) for splitter in settings.REACTION_DATASET_SPLITTERS}
#TODO: set availability of manual splitting up

TOOL_CHOICES = tuple((key, tuple(tool for tool in library.tools)) for key,library in visitorModules.items()) 

class ModelContainer(models.Model):

    """A class for describing a group of statistical models."""

    class Meta:
        app_label = 'DRP'


    library = models.CharField(
        max_length=200, choices=settings.STATS_MODEL_LIBS)
    tool = models.CharField(
        max_length=200, choices=TOOL_CHOICES)
    splitter = models.CharField(
        max_length=200, choices=settings.REACTION_DATASET_SPLITTERS, blank=True, null=True)


    fully_trained = models.ForeignKey("DRP.StatsModel", null=True)

    def __init__(self, predictors, responses, splitter, library, tool, splitter=None, reactions=None, training=None, test=None):
        super(ModelContainer, self).__init__(splitter=splitter, library=library, tool=tool)
        self.reactions = reactions
        self.predictors = predictors
        self.responses = responses
        self.training = training
        self.test = test

    def clean(self):
        if self.tool not visitorModules[self.library].tools:
            raise ValidationError('Selected tool does not exist in selected library', 'wrong_library')
        if self.splitter is None ^ self.reactions is None:
            raise ValidationError('A full set of reactions must be supplied with a splitter', 'argument_mismatch')
        elif (self.test is None ^ self.training is None) ^ self.testLabel is None:
            raise ValidationError('A test set must be supplied with a training set', 'argument_mismatch') 

    @transaction.atomic
    def save(self, *args, **kwargs):
        super(ModelContainer, self).save(*args, **kwargs)
        modelVisitor = visitorModules[self.library].ModelVisitor(self)
        modelVisitor.predictors = predictors
        modelVisitor.responses = responses
        if self.splitter is not None:
            for training, test in splitters[self.splitter].Splitter().split(self.reactions):#wow
                trainingData = DataSet(name=modelVisitor.tag)
                trainingData.save()
                for reaction in training:
                    DataSetRelation.objects.create(reaction=reaction, dataSet = dataSet)
                testingData = 
                modelVisitor.trainingData = trainingData
                modelVisitor.
                testSet = TestSet(model=modelVisitor.stats_model, name= )
        else:
            modelVisitor.trainingData = self.training
            modelVisitor.testData = self.test 

    def build(self, training_reactions, test_reactions, predictors, responses, debug=False):
        """Constructs, trains, and then tests a ML-model using a ModelVisitor
             of the type (library, tool, splitter, etc.) prescribed to this
             ModelContainer."""

        model.setTestingData(test_reactions)

        model._train()
        model._test()

        return model

    def summarize(self):
        """CAUTION: This is a temporary development function. Do not rely on it. """
        """Return a string containing the Confusion Matrices for all stats_models."""
        summaries = "\nK-Fold Validation:\n"
        for model in self.statsmodel_set.all():
            for d in model.predictsDescriptors.all().downcast():
                summaries += "{}\n".format(d.summarize(model))

        return summaries
