#!/usr/bin/env python
"""A test suite for checking that queryset to csv output works."""

import unittest
from decorators import createsUser, joinsLabGroup, createsChemicalClass, createsCompound
from DRP.models import Compound
import csv
from DRPTestCase import DRPTestCase, runTests
from django.conf import settings
loadTests = unittest.TestLoader().loadTestsFromTestCase


@createsUser('Aslan', 'old_magic')
@joinsLabGroup('Aslan', 'Narnia')
@createsChemicalClass('org', 'organic')
@createsCompound('EtOH', 682, 'org', 'Narnia')
@createsCompound('Pyr', 8904, 'org', 'Narnia')
@createsCompound('AcOH', 171, 'org', 'Narnia')
class CsvOutput(DRPTestCase):

    """Tests the output of CSV."""

    def test_regular(self):
        """Test a regular CSV."""
        fn = '/tmp/' + settings.MAIN_SERVER_USER + '_test_csv.csv'
        with open(fn, 'wb') as csvFile:
            Compound.objects.all().toCsv(csvFile)
        with open(fn, 'rb') as csvFile:
            reader = csv.reader(csvFile, delimiter=",", quotechar='"')
            headerRow = reader.next()
            expectedFields = [
                field.name for field in Compound._meta.fields] + ['chemicalClass_1']
            for item in expectedFields:
                self.assertTrue(item in headerRow,
                                'item {} not in {}'.format(item, headerRow))
            self.assertEqual(len(headerRow), len(
                expectedFields), 'headerRow: {0}\ncsvFields: {1}'.format(headerRow, expectedFields))
            rowCount = 0
            for row in reader:
                rowCount += 1
                self.assertEqual(len(row), len(headerRow), 'headerRow ({}): {}\nrow ({}):{}'.format(
                    len(headerRow), headerRow, len(row), row))
            self.assertEqual(rowCount, 3)

    def test_expanded(self):
        """Test an expanded CSV."""
        fn = '/tmp/' + settings.MAIN_SERVER_USER + '_test_csv.csv'
        with open(fn, 'wb') as csvFile:
            Compound.objects.all().toCsv(csvFile, expanded=True)
        with open(fn, 'rb') as csvFile:
            reader = csv.reader(csvFile, delimiter=",")
            headerRow = reader.next()
            expectedFields = [field.name for field in Compound._meta.fields] + [
                'chemicalClass_1'] + [d.csvHeader for d in Compound.objects.all().descriptors]
            for item in expectedFields:
                self.assertTrue(item in headerRow,
                                'item {} not in {}'.format(item, headerRow))
            self.assertEqual(len(headerRow), len(
                expectedFields), 'headerRow: {0}\ncsvFields: {1}'.format(headerRow, expectedFields))
            rowCount = 0
            for row in reader:
                rowCount += 1
                self.assertEqual(len(row), len(headerRow), 'headerRow ({}): {}\nrow ({}):{}'.format(
                    len(headerRow), headerRow, len(row), row))
            self.assertEqual(rowCount, 3)

suite = unittest.TestSuite([
    loadTests(CsvOutput)
])

if __name__ == '__main__':
    runTests(suite)
    # Runs the test- a good way to check that this particular test set works
    # without having to run all the tests.
