#!/usr/bin/env python

#This file contains a (very) loose framework from which others can base their test files and be
#conformant with the local arrangement of test cases.
#For more information about the structure of tests, consult the python documentation at
#https://docs.python.org/2/library/unittest.html

#REMEMBER when you create a new test case to add it to the suite() method, and then
#to have that suite method called in AllTests.py

import unittest
from DRPTestCase import DRPTestCase, runTests
loadTests = unittest.TestLoader().loadTestsFromTestCase

class RidiculousTestCaseOne(DRPTestCase):
    #This class exemplifies the standard structure of a test. Check the documentation for 'rolling your own'

    def setUp(self):
        pass

    def runTest(self):
        pass

    def tearDown(self):
        pass

class RidiculousTestCaseTwo(DRPTestCase):

    def setUp(self):
        pass

    def test_ripeness(self):
        pass

    def tearDown(self):
        pass

suite = unittest.TestSuite([
            loadTests(RidiculousTestCaseOne),
            loadTests(RidiculousTestCaseTwo)
            ])

if __name__=='__main__':
    runTests(suite)
    #Runs the test- a good way to check that this particular test set works without having to run all the tests.
