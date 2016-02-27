from django.conf import settings
import uuid
from DRP.ml_models.model_visitors.weka.AbstractWekaModelVisitor import AbstractWekaModelVisitor
import os

    
class J48(AbstractWekaModelVisitor):

    maxResponseCount = 1

    def __init__(self, *args, **kwargs):
        super(J48, self).__init__(*args, **kwargs)


    def wekaTrainCommand(self, arff_file, filePath, response_index):
        command = "java weka.classifiers.trees.J48 -t {} -d {} -p 0 -c {}".format(arff_file, filePath, response_index)
        return command
        
    def wekaPredictCommand(self, arff_file, model_file, response_index, results_path):
        command = "java weka.classifiers.trees.J48 -T {} -l {} -p 0 -c {} 1> {}".format(arff_file, model_file, response_index, results_path)
        return command
