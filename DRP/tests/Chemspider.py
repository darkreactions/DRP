#!/usr/bin/env python

#REMEMBER when you create a new test case to add it to the suite() method, and then
#to have that suite method called in AllTests.py
import sys
import TestConfig
import DRP.chemspider as chemspider
import DRP.models.CompoundEntry as CompoundEntry
import unittest


loadTests = unittest.TestLoader.loadTestsFromTestCase

class ChemspiderFindFromName(unittest.TestCase):
    #This class exemplifies the standard structure of a test. Check the documentation for 'rolling your own'

    def runTest(self):
        ethanol = chemspider.chemspider_find(['ethanol'])
	self.assertEqual(ethanol.csid,'682')	

class ChemspiderFindFromCAS(unittest.TestCase):

    def runTest(self):
        pyrazole = chemspider.chemspider_find(['288-13-1'])
        self.assertEqual(pyrazole.csid, '1019')

class ChemspiderSearchOnList(unittest.TestCase):

    def runTest(self):
        pyrazole = chemspider.search_chemspider(['288-13-1','pyrazole'])
        self.assertEqual(pyrazole.csid, '1019')

class ChemspiderSearchOnDict(unittest.TestCase):

    def runTest(self):
        ethanol = chemspider.search_chemspider({'compound':'ethanol', 'CAS_ID':'64-17-5'})
        self.assertEqual(ethanol.csid, '682')

class ChemspiderSearchOnEntry(unittest.TestCase):

    def runTest(self):
        c = CompoundEntry.objects.get(pk=1)
        chemspider.search_chemspider(c)
        
class ChemspiderConflict(unittest.TestCase):

    def runTest(self):
        with self.assertRaises(chemspider.ChemspiderError):
            chemspider.search_chemspider(['288-13-1', 'ethanol'])

def suite():
    #This function should be adjusted to contain the loadTests() function enacted on each test case.
    return unittest.TestSuite([
            loadTests(ChemspiderFindFromName),
            loadTests(ChemspiderFindFromCAS),
            loadTests(ChemspiderSearchOnList),
            loadTests(ChemspiderSearchOnDict),
            loadTests(ChemspiderSearchOnEntry),
            loadTests(ChempsiderConflict)
            ])


if __name__ == '__main__':
    #Runs the test- a good way to check that this particular test set works without having to run all the tests.
    unittest.main()
