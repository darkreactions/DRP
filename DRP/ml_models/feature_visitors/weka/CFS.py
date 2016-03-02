from django.conf import settings
import uuid
from DRP.ml_models.model_visitors.weka.AbstractWekaModelVisitor import AbstractWekaModelVisitor
import os


class CFS(AbstractWekaModelVisitor):

    maxResponseCount = 1
    
    def train(self, reactions, descriptorHeaders, filePath):
        arff_file = self._prepareArff(reactions, descriptorHeaders)
        
        # Currently, we support only one "response" variable.
        headers = [h for h in reactions.expandedCsvHeaders if h in descriptorHeaders]
        response_index = headers.index(list(self.statsModel.container.outcomeDescriptors)[0].csvHeader) + 1
        command = "java weka.attributeSelection.CfsSubsetEval -s weka.attributeSelection.BestFirst -i {} -c {}".format(arff_file, filePath, response_index)
        
        feature_selection_output = self._runWekaCommand(command)
    
    def pull_features(self, feature_selection_output):
        
