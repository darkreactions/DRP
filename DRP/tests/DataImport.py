#!/usr/bin/env python
"A test for making sure that data imports from a remote server."

import unittest
from DRPTestCase import DRPTestCase, runTests
from django.db import transaction
from decorators import createsPerformedReaction, createsCompound, joinsLabGroup, createsChemicalClass
from decorators import createsUser, createsCompoundRole, createsRxnDescriptor
from decorators import createsPerformedReactionSet
from decorators import signsExampleLicense 
from DRP.management.commands import import_data
from DRP.models import PerformedReaction
loadTests = unittest.TestLoader().loadTestsFromTestCase

@createsUser('Rorschach', 'whatareyouwaitingfor', is_superuser=True)
@joinsLabGroup('Rorschach', 'Watchmen')
@signsExampleLicense("Rorschach")
@createsPerformedReactionSet
class ApiV1(DRPTestCase):
    """Tests the version 1 api"""

    @transaction.atomic
    def databaseOperation(self, limit=None):
        """Do the actual database import- this is done inside a transaction so that the database does not interfere with itself."""
        c=import_data.Command(limit=None)
        c.handle()

    def test_all(self):
        """Test the database import. For the moment the only test done is to make sure that no exceptions are thrown, since
        any bugs we can identify so far are caught by the database integrity protection from django, which is well tested"""
        self.databaseOperation()
    
    def test_limit(self):
        """test the limited import"""
        self.databaseOperation(5)
        self.assertEqual(PerformedReaction.objects.all().count(), 5)

@createsUser('Rorschach', 'whatareyouwaitingfor')
@joinsLabGroup('Rorschach', 'Watchmen')
@signsExampleLicense("Rorschach")
@createsPerformedReactionSet
class ApiV1FailBadUser(ApiV1):
    pass

suite = unittest.TestSuite([
          loadTests(ApiV1),
          loadTests(ApiV1FailBadUser),
          ])

if __name__=='__main__':
  runTests(suite)
  #Runs the test- a good way to check that this particular test set works without having to run all the tests.
