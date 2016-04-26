#!/usr/bin/env python
'''A module containing tests for the ModelFactory class'''

import unittest
from DRP.models import PerformedReaction, ModelContainer, Descriptor
from decorators import createsPerformedReactionSetOrd, createsPerformedReactionSetBool
from DRPTestCase import DRPTestCase, runTests
from django.conf import settings
loadTests = unittest.TestLoader().loadTestsFromTestCase

class ModelTest(DRPTestCase):
    splitterOptions = None
    visitorOptions = None
    def runTest(self):
        reactions = PerformedReaction.objects.all()
        predictors = Descriptor.objects.filter(heading="testNumber")
        responses = Descriptor.objects.filter(heading="outcome")
        container = ModelContainer.create(self.modelLibrary, self.modelTool, predictors, responses, splitter=self.splitter,
                                          reactions=reactions, splitterOptions=self.splitterOptions, visitorOptions=self.visitorOptions)

        container.build()
        container.save()
        container.full_clean()

        #TODO: We should test the ModelContainer "predict" method here as well.


# TODO XXX more robust testing
# Some of these use bool and some ord
# Should really do tests with more descriptors

@createsPerformedReactionSetOrd
class WekaSVMKFTest(ModelTest):
    modelLibrary = "weka"
    modelTool = "SVM_PUK"
    splitter = "KFoldSplitter"
    
@createsPerformedReactionSetOrd
class WekaSVMMFTest(ModelTest):
    modelLibrary = "weka"
    modelTool = "SVM_PUK"
    splitter = "MutualInfoSplitter"
    
@createsPerformedReactionSetOrd
class WekaJ48KFTest(ModelTest):
    modelLibrary = "weka"
    modelTool = "J48"
    splitter = "KFoldSplitter"
    
@createsPerformedReactionSetBool
class WekaKNNKFTest(ModelTest):
    modelLibrary = "weka"
    modelTool = "KNN"
    splitter = "KFoldSplitter"
    
@createsPerformedReactionSetOrd
class WekaNBKFTest(ModelTest):
    modelLibrary = "weka"
    modelTool = "NaiveBayes"
    splitter = "KFoldSplitter"


suite = unittest.TestSuite([
                    loadTests(WekaSVMKFTest),
                    loadTests(WekaSVMMFTest),
                    loadTests(WekaJ48KFTest),
                    loadTests(WekaKNNKFTest),
                    loadTests(WekaNBKFTest),
                    ])

if __name__=='__main__':
    runTests(suite)
    #Runs the test- a good way to check that this particular test set works without having to run all the tests.
