"""A module containing only the ModelContainer class."""
from django.db import models
from django.conf import settings
from django.db import transaction
from django.core.exceptions import ValidationError
import importlib
visitorModules = {library:importlib.import_module(library) for library in settings.STATS_MODEL_LIBS}
splitters = {splitterName:importlib.import_module(splitter) for splitter in settings.REACTION_DATASET_SPLITTERS}

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
        max_length=200, choices=settings.REACTION_DATASET_SPLITTERS)


    fully_trained = models.ForeignKey("DRP.StatsModel", null=True)

    def __init__(self, reactions, predictors, responses, splitter, library, tool):
        super(ModelContainer, self).__init__(splitter=splitter, library=library, tool=tool)
        self.reactions = reactions
        self.predictors = predictors
        self.responses = responses

    def clean(self):
        if self.tool not visitorModules[self.library].tools:
            raise ValidationError('Selected tool does not exist in selected library', 'wrong_library')

    @transaction.atomic
    def save(self, *args, **kwargs):
        super(ModelContainer, self).save(*args, **kwargs)
        modelVisitor = visitorModules[self.library].ModelVisitor(self)
        for training, test in splitters[self.splitter].Splitter().split(self.reactions):#wow
            modelVisitor.predictors = predictors
            modelVisitor.responses = responses
            modelVisitor.trainingData = training

    def build(self, training_reactions, test_reactions, predictors, responses, debug=False):
        """Constructs, trains, and then tests a ML-model using a ModelVisitor
             of the type (library, tool, splitter, etc.) prescribed to this
             ModelContainer."""

        model.setTrainingData(training_reactions)
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
