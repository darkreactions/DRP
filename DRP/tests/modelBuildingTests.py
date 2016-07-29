#!/usr/bin/env python
"""A module containing tests for the ModelFactory class."""

import unittest
from DRP.models import PerformedReaction, ModelContainer, Descriptor
from .decorators import createsPerformedReactionSetOrd, createsPerformedReactionSetBool
from .drpTestCase import DRPTestCase, runTests
from django.conf import settings
from django.core.management import call_command

loadTests = unittest.TestLoader().loadTestsFromTestCase


class ModelTest(DRPTestCase):
    """Tests for model building."""

    splitterOptions = None
    visitorOptions = None

    def runTest(self):
        """The actual test."""
        reactions = PerformedReaction.objects.all()
        predictors = Descriptor.objects.filter(heading="testNumber")
        responses = Descriptor.objects.filter(heading="outcome")
        container = ModelContainer.create(self.modelLibrary, self.modelTool, predictors, responses, splitter=self.splitter,
                                          reactions=reactions, splitterOptions=self.splitterOptions, visitorOptions=self.visitorOptions)

        container.build()
        container.save()
        container.full_clean()


# Some of these use bool and some ord
# Should really do tests with more descriptors

@createsPerformedReactionSetOrd
class WekaSVMKFTest(ModelTest):
    """Tests Weka SVM."""

    modelLibrary = "weka"
    modelTool = "SVM_PUK"
    splitter = "KFoldSplitter"


@createsPerformedReactionSetOrd
class WekaSVMExpTest(ModelTest):
    """Tests Weka SVM."""

    modelLibrary = "weka"
    modelTool = "SVM_PUK"
    splitter = "ExploratorySplitter"


@createsPerformedReactionSetOrd
class WekaJ48KFTest(ModelTest):
    """Tests Weka j48."""

    modelLibrary = "weka"
    modelTool = "J48"
    splitter = "KFoldSplitter"


@createsPerformedReactionSetBool
class WekaKNNKFTest(ModelTest):
    """Tests Weka KNN."""

    modelLibrary = "weka"
    modelTool = "KNN"
    splitter = "KFoldSplitter"


@createsPerformedReactionSetOrd
class WekaNBKFTest(ModelTest):
    """Some kind of Bayes test..."""

    modelLibrary = "weka"
    modelTool = "NaiveBayes"
    splitter = "KFoldSplitter"



class managementCommandTest():
    call_command("build_model", "-p 'reaction_temperature'")


suite = unittest.TestSuite([
    loadTests(WekaSVMKFTest),
    loadTests(WekaSVMExpTest),
    loadTests(WekaJ48KFTest),
    loadTests(WekaKNNKFTest),
    loadTests(WekaNBKFTest),
])

if __name__ == '__main__':
    runTests(suite)
    # Runs the test- a good way to check that this particular test set works
    # without having to run all the tests.
