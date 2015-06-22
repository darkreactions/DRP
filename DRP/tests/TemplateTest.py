#!/usr/bin/env python

#This file contains a (very) loose framework from which others can base their test files and be
#conformant with the local arrangement of test cases.
#For more information about the structure of tests, consult the python documentation at
#https://docs.python.org/2/library/unittest.html

#REMEMBER when you create a new test case to add it to the suite() method, and then
#to have that suite method called in AllTests.py

import unittest
loadTests = unittest.TestLoader.loadTestsFromTestCase

class RidiculousTestCaseOne(unittest.testcase):

    def setUp(self):
        pass

    def runTest(self):
        pass

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
    return unittest.TestSuite([
            loadTests(RidiculousTestCaseOne),
            loadTests(RidiculousTestCaseTwo)
            ])

if __name__ == '__main__'
    unittest.main()
