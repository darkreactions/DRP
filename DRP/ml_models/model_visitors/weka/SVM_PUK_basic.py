from django.conf import settings
import uuid
from DRP.ml_models.model_visitors.weka.AbstractWekaModelVisitor import AbstractWekaModelVisitor
import os


class SVM_PUK_basic(AbstractWekaModelVisitor):
    wekaCommand = "weka.classifiers.functions.SMO"
    
    def __init__(self, *args, **kwargs):
        super(SVM_PUK_basic, self).__init__(*args, **kwargs)

        self.PUK_OMEGA = 1 #0.5
        self.PUK_SIGMA = 1 #7.0


    #def wekaTrainCommand(self, arff_file, filePath, response_index, weighted=False, cost_matrix_string=None):
        ##kernel = '"weka.classifiers.functions.supportVector.Puk -O {} -S {}"'.format(self.PUK_OMEGA, self.PUK_SIGMA)
        ##if weighted:
            ##if cost_matrix_string is None:
                ##raise RuntimeError("No cost matrix provided but weighted true")
            ##command = "java weka.classifiers.meta.CostSensitiveClassifier -cost-matrix {} -W weka.classifiers.functions.SMO -t {} -d {} -p 0 -c {} -- -K {}".format(cost_matrix_string, arff_file, filePath, response_index, kernel)
        ##else:
            ##command = "java weka.classifiers.functions.SMO -t {} -d {} -K {} -p 0 -c {}".format(arff_file, filePath, kernel, response_index)
        ##return command
        #command = "weka.classifiers.functions.SMO -t {} -d {} -p 0 -c {}".format(arff_file, filePath, response_index)
        #return command

    def wekaTrainOptions(self):
        kernel = '"weka.classifiers.functions.supportVector.Puk -O {} -S {}"'.format(self.PUK_OMEGA, self.PUK_SIGMA)
        return "-K {}".format(kernel)

    #def wekaPredictCommand(self, arff_file, model_file, response_index, results_path):
        #command = "java weka.classifiers.functions.SMO -T {} -l {} -p 0 -c {} 1> {}".format(arff_file, model_file, response_index, results_path)
        #return command
