#!/usr/bin/env python
import sys
import DRP.chemspider as chemspider
import DRP.models.CompoundEntry as CompoundEntry
import unittest
from DRPTestCase import DRPTestCase

loadTests = unittest.TestLoader.loadTestsFromTestCase

class ChemspiderFindFromName(DRPTestCase):
    
    def runTest(self):
        ethanol = chemspider.chemspider_find(['ethanol'])
	self.assertEqual(ethanol.csid,'682')	

class ChemspiderFindFromCAS(DRPTestCase):

    def runTest(self):
        pyrazole = chemspider.chemspider_find(['288-13-1'])
        self.assertEqual(pyrazole.csid, '1019')

class ChemspiderSearchOnList(DRPTestCase):

    def runTest(self):
        pyrazole = chemspider.search_chemspider(['288-13-1','pyrazole'])
        self.assertEqual(pyrazole.csid, '1019')

class ChemspiderSearchOnDict(DRPTestCase):

    def runTest(self):
        ethanol = chemspider.search_chemspider({'compound':'ethanol', 'CAS_ID':'64-17-5'})
        self.assertEqual(ethanol.csid, '682')

class ChemspiderSearchOnEntry(DRPTestCase):

    def runTest(self):
        c = CompoundEntry.objects.get(pk=1)
        chemspider.search_chemspider(c)
        
class ChemspiderConflict(DRPTestCase):

    def runTest(self):
        with self.assertRaises(chemspider.ChemspiderError):
            chemspider.search_chemspider(['288-13-1', 'ethanol'])

def suite():
    return unittest.TestSuite([
            loadTests(ChemspiderFindFromName),
            loadTests(ChemspiderFindFromCAS),
            loadTests(ChemspiderSearchOnList),
            loadTests(ChemspiderSearchOnDict),
            loadTests(ChemspiderSearchOnEntry),
            loadTests(ChempsiderConflict)
            ])


if __name__ == '__main__':
    unittest.main()
