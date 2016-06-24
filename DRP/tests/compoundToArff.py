#!/usr/bin/env python
"""File for testing Arff output."""

import unittest
from .decorators import createsUser, joinsLabGroup, createsChemicalClass
from .decorators import createsCompound
from DRP.models import Compound
from django.conf import settings
import subprocess
from .drpTestCase import DRPTestCase, runTests

loadTests = unittest.TestLoader().loadTestsFromTestCase


@createsUser('Aslan', 'old_magic')
@joinsLabGroup('Aslan', 'Narnia')
@createsChemicalClass('Org', 'organic')
@createsCompound('2-amep', 104820, 'Org', 'Narnia', custom=False)
@createsCompound('hmta', 3959, 'Org', 'Narnia', custom=False)
@createsCompound('sp', 1071, 'Org', 'Narnia', custom=False)
class CompoundToArff(DRPTestCase):

    """Validates the structure of the Arff- this could be more detailed if we find that we encounter issues later."""

    def checkArff(self, fn):
        """Check the arff is valid."""
        process = subprocess.Popen(['java', '-cp', settings.WEKA_PATH['3.6'],
                                    'weka.core.Instances', fn], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        process.wait()
        c = process.returncode
        # on the off chance weka ever returns a non-zero error code
        self.assertEqual(0, c)
        res, resErr = process.communicate()
        self.assertFalse(resErr)

    def test_regular(self):
        """Test a regular arff."""
        fn = '/tmp/' + settings.MAIN_SERVER_USER + '_test_csv.arff'
        with open(fn, 'wb') as arffFile:
            Compound.objects.all().toArff(arffFile)
        self.checkArff(fn)

    def test_expanded(self):
        """Test an expanded arff."""
        fn = '/tmp/' + settings.MAIN_SERVER_USER + '_test_ex_csv.arff'
        with open(fn, 'wb') as arffFile:
            Compound.objects.all().toArff(arffFile, expanded=True)
        self.checkArff(fn)

suite = unittest.TestSuite([
    loadTests(CompoundToArff)
])

if __name__ == '__main__':
    runTests(suite)
    # Runs the test- a good way to check that this particular test set works
    # without having to run all the tests.
