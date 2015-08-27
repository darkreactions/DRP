#!/usr/bin/env python
'''A test suite for checking that queryset to csv output works'''

import unittest
from DRPTestCase import DRPTestCase, runTests
from decorators import createsUser, joinsLabGroup, createsChemicalClass, createsCompound
from DRP.models import Compound, MolDescriptor
import csv
loadTests = unittest.TestLoader().loadTestsFromTestCase

@createsUser('Aslan', 'old_magic')
@joinsLabGroup('Aslan', 'Narnia')
@createsChemicalClass('org', 'organic')
@createsCompound('EtOH', 682, 'org', 'Narnia')
@createsCompound('Pyr', 8904, 'org', 'Narnia')
@createsCompound('AcOH', 171, 'org', 'Narnia')
class CsvOutput(DRPTestCase):
  #This class exemplifies the standard structure of a test. Check the documentation for 'rolling your own'

  def test_regular(self):
    fn = '/tmp/test_csv.csv'
    with open(fn, 'wb') as csvFile:
      Compound.objects.all().toCsv(csvFile)
    with open(fn, 'rb') as csvFile:
      reader = csv.reader(csvFile, delimiter=",", quotechar='"')
      headerRow = reader.next()
      expectedFields=[field.name for field in Compound._meta.fields] + ['chemicalClass_1']
      for item in expectedFields:
        self.assertTrue(item in headerRow, 'item {} not in {}'.format(item, headerRow))
      self.assertEqual(len(headerRow), len(expectedFields), 'headerRow: {0}\ncsvFields: {1}'.format(headerRow, expectedFields))
      rowCount = 0
      for row in reader:
        rowCount +=1
        self.assertEqual(len(row), len(headerRow), 'headerRow ({}): {}\nrow ({}):{}'.format(len(headerRow),headerRow,len(row), row))
      self.assertEqual(rowCount, 3)

  def test_expanded(self):
    fn = '/tmp/test_ex_csv.csv'
    with open(fn, 'wb') as csvFile:
      Compound.objects.all().toCsv(csvFile, expanded=True) 
    with open(fn, 'rb') as csvFile:
      reader = csv.reader(csvFile, delimiter=",")
      headerRow = reader.next()
      expectedFields=[field.name for field in Compound._meta.fields] + ['chemicalClass_1'] + [d.csvHeader for d in Compound.objects.all().descriptors()]
      for item in expectedFields:
        self.assertTrue(item in headerRow, 'item {} not in {}'.format(item, headerRow))
      self.assertEqual(len(headerRow), len(expectedFields), 'headerRow: {0}\ncsvFields: {1}'.format(headerRow, expectedFields))
      rowCount = 0
      for row in reader:
        rowCount += 1
        self.assertEqual(len(row), len(headerRow), 'headerRow ({}): {}\nrow ({}):{}'.format(len(headerRow),headerRow,len(row), row))
      self.assertEqual(rowCount, 3)
  

suite = unittest.TestSuite([
          loadTests(CsvOutput)
          ])

if __name__=='__main__':
  runTests(suite)
  #Runs the test- a good way to check that this particular test set works without having to run all the tests.
