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

        mod = self.getVisitorModule()
        model = mod.ModelVisitor(self)

        model.setPredictors(predictors)
        model.setResponses(responses)
        model.setTrainingData(training_reactions)
        model.setTestingData(test_reactions)

        model._train()
        model._test()

        return model

    def summarize(self):
        """CAUTION: This is a temporary development function. Do not rely on it. """
        """Return a string containing the Confusion Matrices for all stats_models."""
        summaries = "\nK-Fold Validation:\n"
        mod = self.getVisitorModule()
        for model in self.statsmodel_set.all():
            preds = mod.ModelVisitor(self, stats_model=model).getConfusionMatrices()
            summaries += "{}\n".format(preds)

        return summaries

    def getVisitorModule(self):
        """Return the module of the ModelVisitor that can wrap the stats_models
           owned by this ModelContainer."""
        module_name = "{}_{}".format(self.library, self.tool)
        try:
            module_path = "DRP.ml_models.model_visitors.{}".format(module_name)
            module = importlib.import_module(module_path)
        except ImportError:
            error = "Cannot import model visitor: \"{}\".".format(module_name)
            raise NotImplementedError(error)
        return module
