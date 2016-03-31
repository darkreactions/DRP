from django.conf import settings
import uuid
from DRP.ml_models.model_visitors.weka.AbstractWekaModelVisitor import AbstractWekaModelVisitor
import os


class SVM_PUK_basic(AbstractWekaModelVisitor):
    def __init__(self, *args, **kwargs):
        super(SVM_PUK_basic, self).__init__(*args, **kwargs)

        self.PUK_OMEGA = 1 #0.5
        self.PUK_SIGMA = 1 #7.0

    def wekaTrainCommand(self, arff_file, filePath, response_index):
        kernel = "\"weka.classifiers.functions.supportVector.Puk -O {} -S {}\"".format(self.PUK_OMEGA, self.PUK_SIGMA)
        command = "java weka.classifiers.functions.SMO -t {} -d {} -K {} -p 0 -c {}".format(arff_file, filePath, kernel, response_index)
        return command

    def wekaPredictCommand(self, arff_file, model_file, response_index, results_path):
        command = "java weka.classifiers.functions.SMO -T {} -l {} -p 0 -c {} 1> {}".format(arff_file, model_file, response_index, results_path)
        return command
