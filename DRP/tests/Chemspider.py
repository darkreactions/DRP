#!/usr/bin/env python

#REMEMBER when you create a new test case to add it to the suite() method, and then
#to have that suite method called in AllTests.py
import TestConfig
import chemspider
import unittest
loadTests = unittest.TestLoader.loadTestsFromTestCase

class ChemspiderFind (unittest.testcase):
    #This class exemplifies the standard structure of a test. Check the documentation for 'rolling your own'

    def setUp(self):
        pass

    def runTest(self):
        chemspider.chemspider_find(['ethanol'])
	

    def tearDown(self):
        pass

class RidiculousTestCaseTwo(unittest.testcase):

    def setUp(self):
        pass

    def test_ripeness(self):
        pass

    def tearDown(self):
        pass

def suite():
    #This function should be adjusted to contain the loadTests() function enacted on each test case.
    return unittest.TestSuite([
            loadTests(RidiculousTestCaseOne),
            loadTests(RidiculousTestCaseTwo)
            ])

if __name__ == '__main__'
    #Runs the test- a good way to check that this particular test set works without having to run all the tests.
    unittest.main()
