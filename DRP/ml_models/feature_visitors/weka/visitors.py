from AbstractWekaFeatureVisitor import AbstractWekaFeatureVisitor
import importlib
from django.conf import settings

modelVisitorModules = {library:importlib.import_module(settings.STATS_MODEL_LIBS_DIR + "."+ library) for library in settings.STATS_MODEL_LIBS}

class CFS(AbstractWekaFeatureVisitor):
    wekaCommand = "weka.attributeSelection.CfsSubsetEval -s weka.attributeSelection.BestFirst"

class InfoGain(AbstractWekaFeatureVisitor):
    wekaCommand = "weka.attributeSelection.InfoGainAttributeEval -s weka.attributeSelection.Ranker"
    
class ChiSquared(AbstractWekaFeatureVisitor):
    wekaCommand = "weka.attributeSelection.ChiSquaredAttributeEval -s weka.attributeSelection.Ranker"

class Wrapper(AbstractWekaFeatureVisitor):
    def __init__(self, modelVisitorTool="J48", modelVisitorLibrary="weka", modelVisitorOptions={}, *args, **kwargs):
        getattr(modelVisitorModules[modelVisitorLibrary], modelVisitorTool)
        modelVisitor = getattr(modelVisitorModules[modelVisitorLibrary], modelVisitorTool)(statsModel=None, **modelVisitorOptions)
        self.wekaCommand = "weka.attributeSelection.WrapperSubsetEval -s weka.attributeSelection.BestFirst -B {}".format(modelVisitor.wekaCommand)
        self.wekaOptions = modelVisitor.wekaTrainOptions
        
        super(Wrapper, self).__init__(*args, **kwargs)

