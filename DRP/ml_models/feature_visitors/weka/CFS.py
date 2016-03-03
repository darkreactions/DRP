from django.conf import settings
import uuid
from AbstractWekaFeatureVisitor import AbstractWekaFeatureVisitor
import os


class CFS(AbstractWekaFeatureVisitor):

    maxResponseCount = 1
    
    def wekaTrainCommand(self, arff_file, response_index):
        command = "java weka.attributeSelection.CfsSubsetEval -s weka.attributeSelection.BestFirst -i {} -c {}".format(arff_file, response_index)
        return command
