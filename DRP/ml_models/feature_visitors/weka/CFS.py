from django.conf import settings
import uuid
from DRP.ml_models.model_visitors.weka.AbstractWekaModelVisitor import AbstractWekaModelVisitor
import os


class CFS(AbstractWekaModelVisitor):

    maxResponseCount = 1
    
    def wekaTrainCommand(self, arff_file, response_index):
        command = "java weka.attributeSelection.CfsSubsetEval -s weka.attributeSelection.BestFirst -i {} -c {} &> {}".format(arff_file, response_index)
