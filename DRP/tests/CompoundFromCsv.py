#!/usr/bin/env python

"""Tests for ingesting CSVs into the  compound guide."""

# This file contains a (very) loose framework from which others can base their test files and be
# conformant with the local arrangement of test cases.
# For more information about the structure of tests, consult the python documentation at
# https://docs.python.org/2/library/unittest.html

# REMEMBER when you create a new test case to add it to the suite() method, and then
# to have that suite method called in AllTests.py

import unittest
from DRPTestCase import DRPTestCase, runTests
from DRP.models import Compound, LabGroup
from decorators import createsUser, joinsLabGroup, createsChemicalClass
import os.path
from django.core.exceptions import ValidationError
from django.conf import settings
loadTests = unittest.TestLoader().loadTestsFromTestCase


@createsUser('Aslan', 'testingpass')
@joinsLabGroup('Aslan', 'Narnia')
@createsChemicalClass('Org', 'Organic')
class Good(DRPTestCase):
    """Tests the spreadsheets that should work (whose names end with the ssids)."""

    ssids = (1, 3, 5, 9, 11, 13)
    filenameStub = 'compound_spread_test{0}.csv'
    prefix = os.path.join(settings.APP_DIR, os.path.join('tests', 'resource'))

    @property
    def fileNames(self):
        for id in self.ssids:
            yield os.path.join(self.prefix, self.filenameStub.format(id))

    def runTest(self):
        for filename in self.fileNames:
            compounds = Compound.objects.fromCsv(filename, LabGroup.objects.get(title='Narnia'))
            for compound in compounds:
                compound.csConsistencyCheck()
                compound.full_clean()
                compound.save()
            self.assertEqual(len(compounds), 8)
            for compound in compounds:
                compound.delete()


class Broken(Good):
    """Tests the broken spreads whose names ends with the values in ssids."""

    ssids = (2, 4, 6, 7, 8, 10, 12, 14, 15)

    def runTest(self):
        for fileName in self.fileNames:
            with self.assertRaises(ValidationError):
                compounds = Compound.objects.fromCsv(fileName, LabGroup.objects.get(title='Narnia'))
                for compound in compounds:
                    compound.csConsistencyCheck()
                    compound.full_clean()
                    compound.save()
                Compound.objects.all().delete()

suite = unittest.TestSuite([
    loadTests(Good),
    loadTests(Broken)
])

if __name__ == '__main__':
    runTests(suite)
    # Runs the test- a good way to check that this particular test set works without having to run all the tests.
