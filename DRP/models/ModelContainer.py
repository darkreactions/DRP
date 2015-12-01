"""A module containing only the ModelContainer class."""
from django.db import models
from django.conf import settings
import importlib

class ModelContainer(models.Model):

    """A class for describing a group of statistical models."""

    class Meta:
        app_label = 'DRP'


    library = models.CharField(
        max_length=200, choices=settings.LIBRARY_CHOICES)
    tool = models.CharField(
        max_length=200, choices=settings.TOOL_CHOICES)
    splitter = models.CharField(
        max_length=200, choices=settings.SPLITTER_CHOICES)


    fully_trained = models.ForeignKey("DRP.StatsModel", null=True)

    def build(self, training_reactions, test_reactions, predictors, responses, debug=True):
        """Constructs, trains, and then tests a ML-model using a ModelVisitor
             of the type (library, tool, splitter, etc.) prescribed to this
             ModelContainer."""

        if debug: print "PREPPING - 1"
        mod = self.getVisitorModule()
        model = mod.ModelVisitor(self)

        if debug: print "PREPPING - 2"

        model.setPredictors(predictors)
        model.setResponses(responses)

        if debug: print "PREPPING - 3"
        model.setTrainingData(training_reactions)
        model.setTestingData(test_reactions)

        if debug: print "TRAINING"
        model._train()
        if debug: print "TESTING"
        model._test()

        return model

    def summarize(self):
        """CAUTION: This is a temporary development function. Do not rely on it. """
        """Return a string containing the Confusion Matrices for all stats_models."""
        summaries = "\nK-Fold Validation:\n"
        for model in self.statsmodel_set.all():
            for d in model.predictsDescriptors.all().downcast():
                summaries += d.summarize(model) + "\n"

        return summaries

    def getVisitorModule(self):
        """Return the module of the ModelVisitor that can wrap the stats_models
           owned by this ModelContainer."""
        module_name = "{}_{}".format(self.library, self.tool)
        try:
            module_path = "DRP.ml_models.model_visitors.{}".format(module_name)
            module = importlib.import_module(module_path)
        except ImportError as e:
            error = "Cannot import model visitor: \"{}\": {}.".format(module_name, e)
            raise NotImplementedError(error)
        return module
