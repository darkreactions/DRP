#!/usr/bin/env python
'''A module containing tests for the ModelFactory class'''

import unittest
from DRP.models import PerformedReaction
from DRP.ml_models.ModelFactory import ModelFactory
from decorators import createsPerformedReaction, createsCompound, joinsLabGroup, createsChemicalClass, createsUser, createsCompoundRole, createsRxnDescriptor
from DRPTestCase import DRPTestCase, runTests
loadTests = unittest.TestLoader().loadTestsFromTestCase

@createsUser('Rorschach', 'whatareyouwaitingfor')
@joinsLabGroup('Rorschach', 'Watchmen')
@createsChemicalClass('Org', 'Organic')
@createsChemicalClass('Water', 'Water')
@createsCompound('EtOH', 682, 'Org', 'Watchmen')
@createsCompound('Water', 937, 'Water', 'Watchmen')
@createsCompoundRole('Org', 'Organic')
@createsCompoundRole('Water', 'Water')
@createsRxnDescriptor("outcome", "OrdRxnDescriptor", {"maximum":4, "minimum":1})
@createsPerformedReaction("Watchmen", "Rorschach",["EtOH"], ["Org"], [0.1], {"outcome":1})
@createsPerformedReaction("Watchmen", "Rorschach",["EtOH"], ["Org"], [0.2],{"outcome":1})
@createsPerformedReaction("Watchmen", "Rorschach",["EtOH"], ["Org"], [0.3],{"outcome":1})
@createsPerformedReaction("Watchmen", "Rorschach",["EtOH", "Water"], ["Org", "Water"],[0.3, 0.3],{"outcome":3})
@createsPerformedReaction("Watchmen", "Rorschach",["EtOH", "Water"], ["Org", "Water"],[0.4, 0.5],{"outcome":4})
@createsPerformedReaction("Watchmen", "Rorschach",["Water", "EtOH"], ["Water", "Org"],[0.4, 0.5],{"outcome":4})
class BasicWekaSVM(DRPTestCase):

  def runTest(self):
    reactions = PerformedReaction.objects.all()

    #TODO: Test Headers
    headers = ["compound_0_amount", "compound_1_amount", "valid"]

    factory = ModelFactory()
    factory.build(reactions, headers, modelLibrary="weka", modelType="svm",
                  debug=True)


suite = unittest.TestSuite([
          loadTests(BasicWekaSVM),
          ])

if __name__=='__main__':
  runTests(suite)
  #Runs the test- a good way to check that this particular test set works without having to run all the tests.
