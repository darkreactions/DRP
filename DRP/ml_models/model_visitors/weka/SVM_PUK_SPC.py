from django.conf import settings
import uuid
from DRP.ml_models.model_visitors.weka.AbstractWekaModelVisitor import AbstractWekaModelVisitor
from DRP.models.descriptors import BooleanDescriptor, NumericDescriptor, CategoricalDescriptor, OrdinalDescriptor
from DRP.models.rxnDescriptorValues import BoolRxnDescriptorValue, OrdRxnDescriptorValue, BoolRxnDescriptorValue
import os


class SVM_PUK_SPC(AbstractWekaModelVisitor):
    """
    Optimizes for specificity.
    Therefore completely ignores true data points.
    Only works for two classes. Useful for testing only.
    """

    maxResponseCount = 1

    def __init__(self, *args, **kwargs):
        super(SVM_PUK_SPC, self).__init__(*args, **kwargs)

        self.PUK_OMEGA = 0.5
        self.PUK_SIGMA = 7.0


    def wekaTrainCommand(self, arff_file, filePath, response_index):
        raise RuntimeError("Because the SVM optimized for SPC must know the class counts, "
                            "it implements its own train function and does not call wekaTrain.")
                            
    def train(self, reactions, descriptorHeaders, filePath, verbose=False):
        arff_file = self._prepareArff(reactions, descriptorHeaders)

        ## Currently, we support only one "response" variable.
            
        ## the i,j entry in cost matrix corresponds to cost of classifying an instance of class i as class j
        ## since we want misclassification of an instance to be equally costly regardless of what it is classified as
        ## all entries in a row not on the diagonal should be the same. The diagonal should always be 0 as correct
        ## classification has no cost. So that every class is weighted the same as a whole, each class's weight is 
        ## 1/(number of instances of that class). To reduce floating point arithmetic errors, we first multiply by
        ## total number of data points, so each class is weighted by (total_instances)/(class_count)
        ## classes for which class_count is 0 do not matter, so their weight is 0 (to avoid division by 0)
        ## For boolean classification, True is class 0 and False is class 1 (Because that's how it's set up in the toArff function)
        ## TODO XXX Actually check the order of classes from the arff file?

        cost_matrix_string = '"[0.0 0.0; 1.0 0.0]"'
     
        # Currently, we support only one "response" variable.
        headers = [h for h in reactions.expandedCsvHeaders if h in descriptorHeaders]
        response_index = headers.index(list(self.statsModel.container.outcomeDescriptors)[0].csvHeader) + 1
        
        kernel = '"weka.classifiers.functions.supportVector.Puk -O {} -S {}"'.format(self.PUK_OMEGA, self.PUK_SIGMA)
        command = "java weka.classifiers.meta.CostSensitiveClassifier -cost-matrix {} -W weka.classifiers.functions.SMO -t {} -d {} -p 0 -c {} -- -K {}".format(cost_matrix_string, arff_file, filePath, response_index, kernel)
        self._runWekaCommand(command)

    def wekaPredictCommand(self, arff_file, model_file, response_index, results_path):
        command = "java weka.classifiers.functions.SMO -T {} -l {} -p 0 -c {} 1> {}".format(arff_file, model_file, response_index, results_path)
        return command

