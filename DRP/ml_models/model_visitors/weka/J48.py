from django.conf import settings
import uuid
from DRP.ml_models.model_visitors.weka.AbstractWekaModelVisitor import AbstractWekaModelVisitor
import os

    
class J48(AbstractWekaModelVisitor):
    wekaCommand = "weka.classifiers.trees.J48"
