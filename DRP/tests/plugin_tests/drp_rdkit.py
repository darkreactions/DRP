"""A test suite for the drp_rdkit plugin for calculating molecular descriptors."""

import unittest
import math
import sys
from ..drpTestCase import DRPTestCase, runTests
from DRP.tests.decorators import createsUser, joinsLabGroup, createsChemicalClass
from DRP.tests.decorators import createsCompound
#from DRP.plugins.moldescriptors import drp_rdkit
from DRP.tests.decorators import createsCompound
from DRP.models import Compound, NumMolDescriptorValue, BoolMolDescriptorValue
import logging
tracer = logging.getLogger('DRP.tracer')
tracer.debug('Importing Complete')
loadTests = unittest.TestLoader().loadTestsFromTestCase


@createsUser('Aslan', 'old_magic')
@joinsLabGroup('Aslan', 'Narnia')
@createsChemicalClass('inorg', 'inorganic')
@createsCompound('VOx', 14130, 'inorg', 'Narnia')
class DescriptorCalculation(DRPTestCase):
    """Tests for descriptor calculations. Presently rudimentary."""

    deleteDescriptors = False

    def test_molecularWeight(self):
        """Test for the molecular weight calculation."""
        tracer.debug('MWTest')
        nums, bools = drp_rdkit.calculate(
            Compound.objects.get(CSID=14130), whitelist=['mw'])
        self.assertTrue(math.fabs(
            nums[0].value - 181.8800) < 0.001, math.fabs(nums[0].value - 181.8800))
        self.assertEqual(len(bools), 0)
        drp_rdkit.descriptorDict.initialised = False

    def test_rbc(self):
        """Test for the rotatable bond count calculation."""
        tracer.debug('RBCTest')
        nums, bools = drp_rdkit.calculate(
            Compound.objects.get(CSID=14130), whitelist=['rbc'])
        self.assertEqual(len(nums), 1)
        self.assertEqual(nums[0].descriptor.heading, 'rbc')
        drp_rdkit.descriptorDict.initialised = False

    def test_chi0v(self):
        """Test for teh Chi0v calculation."""
        tracer.debug('ChiTest')
        nums, bools = drp_rdkit.calculate(
            Compound.objects.get(CSID=14130), whitelist=['Chi0v'])
        self.assertEqual(len(nums), 1)
        self.assertEqual(nums[0].descriptor.heading, 'Chi0v')
        drp_rdkit.descriptorDict.initialised = False

    def test_ox(self):
        """Test the boolean descriptors for the presence of atoms in certain oxidation states."""
        tracer.debug('OxTest')
        nums, bools = drp_rdkit.calculate(
            Compound.objects.get(CSID=14130), whitelist=['V@5'])
        self.assertTrue(bools[0].value)
        self.assertEqual(len(nums), 0)
        drp_rdkit.descriptorDict.initialised = False

    def test_all(self):
        """Test the calculation of all descriptors."""
        tracer.debug('AllTest')
        nums, bools = drp_rdkit.calculate(Compound.objects.get(CSID=14130))
        self.assertEqual(len(nums), 3)
        self.assertEqual(len(bools), 405)
        drp_rdkit.descriptorDict.initialised = False


@createsUser('Aslan', 'old_magic')
@joinsLabGroup('Aslan', 'Narnia')
@createsChemicalClass('org', 'Organic')
@createsChemicalClass('inorg', 'inorganic')
@createsCompound('VOx', 14130, 'inorg', 'Narnia')
@createsCompound('EtOH', 682, 'org', 'Narnia')
class MultipleDescriptorCalculation(DRPTestCase):
    """Tests for the calculation of multiple descriptors simultaneously."""

    deleteDescriptors = False

    def test_all(self):
        """Test the calculation of all descriptors."""
        tracer.debug('AllTest2')
        drp_rdkit.calculate_many(Compound.objects.all())
        tracer.debug("Here is a list: {}".format(
            str(NumMolDescriptorValue.objects.all())))
        self.assertEqual(NumMolDescriptorValue.objects.all().count(), 6)
        self.assertEqual(BoolMolDescriptorValue.objects.all().count(), 810, ', '.join(
            str(dvo) + '\n' for dvo in BoolMolDescriptorValue.objects.all()))
        drp_rdkit.descriptorDict.initialised = False

    def test_whitelist(self):
        """Test the calculation of a whitelist of descriptors."""
        tracer.debug('WLTest')
        drp_rdkit.calculate_many(Compound.objects.all(), whitelist=['mw'])
        self.assertEqual(NumMolDescriptorValue.objects.all().count(), 2)
        self.assertEqual(BoolMolDescriptorValue.objects.all().count(), 0)
        drp_rdkit.descriptorDict.initialised = False

suite = unittest.TestSuite([
    loadTests(DescriptorCalculation),
    loadTests(MultipleDescriptorCalculation)
])
