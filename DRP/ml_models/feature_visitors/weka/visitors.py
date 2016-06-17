"""Feature visitors appropriate to the weka modelling engine."""
from AbstractWekaFeatureVisitor import AbstractWekaFeatureVisitor
import importlib
from django.conf import settings

modelVisitorModules = {library: importlib.import_module(settings.STATS_MODEL_LIBS_DIR + "." + library) for library in settings.STATS_MODEL_LIBS}


class CFS(AbstractWekaFeatureVisitor):

    """Best first CFS subset evaluation by Weka."""

    wekaCommand = "weka.attributeSelection.CfsSubsetEval -s weka.attributeSelection.BestFirst"


class InfoGain(AbstractWekaFeatureVisitor):

    """Ranked feature selection by information gain."""

    wekaCommand = "weka.attributeSelection.InfoGainAttributeEval -s weka.attributeSelection.Ranker"


class ChiSquared(AbstractWekaFeatureVisitor):

    """Feature selection based on chi-squared distribution evaluations."""

    wekaCommand = "weka.attributeSelection.ChiSquaredAttributeEval -s weka.attributeSelection.Ranker"


class Wrapper(AbstractWekaFeatureVisitor):
    """A wrapper for performing feature visitation with respect to a specific model class."""

    def __init__(self, modelVisitorTool="J48", modelVisitorLibrary="weka", modelVisitorOptions={}, *args, **kwargs):
        """Specify the model class to be used."""
        getattr(modelVisitorModules[modelVisitorLibrary], modelVisitorTool)
        modelVisitor = getattr(modelVisitorModules[modelVisitorLibrary], modelVisitorTool)(statsModel=None, **modelVisitorOptions)
        self.wekaCommand = "weka.attributeSelection.WrapperSubsetEval -s weka.attributeSelection.BestFirst -B {}".format(modelVisitor.wekaCommand)
        self.wekaOptions = modelVisitor.wekaTrainOptions

        super(Wrapper, self).__init__(*args, **kwargs)
