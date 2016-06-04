#!/usr/bin/env python

# This file contains a (very) loose framework from which others can base their test files and be
# conformant with the local arrangement of test cases.
# For more information about the structure of tests, consult the python documentation at
# https://docs.python.org/2/library/unittest.html

# REMEMBER when you create a new test case to add it to the suite() method, and then
# to have that suite method called in AllTests.py

import unittest
from decorators import createsUser, joinsLabGroup
from decorators import loadsCompoundsFromCsv
from DRP.models import Compound
from django.conf import settings
import subprocess
from DRPTestCase import DRPTestCase, runTests

loadTests = unittest.TestLoader().loadTestsFromTestCase


@createsUser('Aslan', 'old_magic')
@joinsLabGroup('Aslan', 'Narnia')
@loadsCompoundsFromCsv('Narnia', 'compound_spread_test1.csv')
class CompoundToArff(DRPTestCase):
    '''Validates the structure of the Arff- this could be more detailed if we find that we encounter issues later'''

    def checkArff(self, fn):
        process = subprocess.Popen(['java', '-cp', settings.WEKA_PATH['3.6'], 'weka.core.Instances', fn], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        process.wait()
        c = process.returncode
        self.assertEqual(0, c)  # on the off chance weka ever returns a non-zero error code
        res, resErr = process.communicate()
        self.assertFalse(resErr)

    def test_regular(self):
        fn = '/tmp/test_csv.arff'
        with open(fn, 'wb') as arffFile:
            Compound.objects.all().toArff(arffFile)
        self.checkArff(fn)

    def test_expanded(self):
        fn = '/tmp/test_ex_csv.arff'
        with open(fn, 'wb') as arffFile:
            Compound.objects.all().toArff(arffFile, expanded=True)
        self.checkArff(fn)

suite = unittest.TestSuite([
    loadTests(CompoundToArff)
])

if __name__ == '__main__':
    runTests(suite)
    # Runs the test- a good way to check that this particular test set works without having to run all the tests.
