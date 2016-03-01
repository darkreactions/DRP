#!/usr/bin/env python
'''A module containing tests for the ModelFactory class'''

import unittest
from DRP.models import PerformedReaction, ModelContainer, Descriptor
from decorators import createsPerformedReactionSetOrd, createsPerformedReactionSetBool
from DRPTestCase import DRPTestCase, runTests
loadTests = unittest.TestLoader().loadTestsFromTestCase

class ModelTest(DRPTestCase):
    featureLibrary = None
    featureTool = None
    
    def runTest(self):
        reactions = PerformedReaction.objects.all()
        container = ModelContainer.create(self.modelLibrary, self.modelTool, splitter=self.splitter,
                                          reactions=reactions)
        container.save()
        predictors = Descriptor.objects.filter(heading="testNumber")
        responses = Descriptor.objects.filter(heading="outcome")
        container.build(predictors, responses)
        container.save()

        #TODO: We should test the ModelContainer "predict" method here as well.


# TODO XXX more robust testing
# Some of these use bool and some ord
# Should really do tests with more descriptors

@createsPerformedReactionSetOrd
class WekaSVMBasicKFTest(ModelTest):
    modelLibrary = "weka"
    modelTool = "SVM_PUK_basic"
    splitter = "KFoldSplitter"
    
@createsPerformedReactionSetOrd
class WekaSVMBasicMFTest(ModelTest):
    modelLibrary = "weka"
    modelTool = "SVM_PUK_basic"
    splitter = "MutualInfoSplitter"
    
@createsPerformedReactionSetBool
class WekaSVMBCRKFTest(ModelTest):
    modelLibrary = "weka"
    modelTool = "SVM_PUK_BCR"
    splitter = "KFoldSplitter"

    
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
                    loadTests(WekaSVMBasicKFTest),
                    loadTests(WekaSVMBasicMFTest),
                    loadTests(WekaSVMBCRKFTest),
                    loadTests(WekaJ48KFTest),
                    loadTests(WekaKNNKFTest),
                    loadTests(WekaNBKFTest),
                    ])

if __name__=='__main__':
    runTests(suite)
    #Runs the test- a good way to check that this particular test set works without having to run all the tests.
