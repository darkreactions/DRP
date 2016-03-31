from django.conf import settings
import uuid
from DRP.ml_models.model_visitors.weka.AbstractWekaModelVisitor import AbstractWekaModelVisitor
from DRP.models.descriptors import BooleanDescriptor, NumericDescriptor, CategoricalDescriptor, OrdinalDescriptor
from DRP.models.rxnDescriptorValues import BoolRxnDescriptorValue, OrdRxnDescriptorValue, BoolRxnDescriptorValue
import os


class SVM_PUK_BCR(AbstractWekaModelVisitor):

    maxResponseCount = 1

    def __init__(self, *args, **kwargs):
        super(SVM_PUK_BCR, self).__init__(*args, **kwargs)

        self.PUK_OMEGA = 1 #0.5
        self.PUK_SIGMA = 1 #7.0


    def wekaTrainCommand(self, arff_file, filePath, response_index):
        # TODO XXX parse the arff to get these instead of throwing an error.
        raise RuntimeError("Because the SVM optimized for BCR must know the class counts, "
                            "it implements its own train function and does not call wekaTrain.")

    def train(self, reactions, descriptorHeaders, filePath, verbose=False, BCR=True):
        arff_file = self._prepareArff(reactions, descriptorHeaders, verbose=verbose)

        ## Currently, we support only one "response" variable.
        response = list(self.statsModel.container.outcomeDescriptors)[0]

        cost_matrix_string = self.BCR_cost_matrix(reactions, response)
     
        # Currently, we support only one "response" variable.
        headers = [h for h in reactions.expandedCsvHeaders() if h in descriptorHeaders]
        response_index = headers.index(response.csvHeader) + 1

        kernel = '"weka.classifiers.functions.supportVector.Puk -O {} -S {}"'.format(self.PUK_OMEGA, self.PUK_SIGMA)
        command = "java weka.classifiers.meta.CostSensitiveClassifier -cost-matrix {} -W weka.classifiers.functions.SMO -t {} -d {} -p 0 -c {} -- -K {}".format(cost_matrix_string, arff_file, filePath, response_index, kernel)

        self._runWekaCommand(command, verbose=verbose)

    def wekaPredictCommand(self, arff_file, model_file, response_index, results_path):
        command = "java weka.classifiers.functions.SMO -T {} -l {} -p 0 -c {} 1> {}".format(arff_file, model_file, response_index, results_path)
        return command

