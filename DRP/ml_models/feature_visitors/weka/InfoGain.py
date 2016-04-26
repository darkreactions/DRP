from django.conf import settings
import uuid
from DRP.ml_models.model_visitors.weka.AbstractWekaModelVisitor import AbstractWekaModelVisitor
import os


class InfoGain(AbstractWekaModelVisitor):

  maxResponseCount = 1

  def __init__(self, *args, **kwargs):
    super(InfoGain, self).__init__(*args, **kwargs)


  def train(self, reactions, descriptorHeaders, filePath):
    arff_file = self._prepareArff(reactions, descriptorHeaders)

    # Currently, we support only one "response" variable.
    headers = [h for h in reactions.expandedCsvHeaders if h in descriptorHeaders]
    response_index = headers.index(list(self.statsModel.container.outcomeDescriptors)[0].csvHeader) + 1
    command = "java weka.attributeSelection.InfoGainAttributeEval -s weka.attributeSelection.Ranker -i {} -c {}".format(arff_file, filePath, response_index)

    feature_selection_output = self._runWekaCommand(command)

  def pull_features(self, feature_selection_output):
    pass
