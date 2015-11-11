#!/usr/bin/env python
"""A module containing tests for molecular descriptor classes."""

import unittest
from DRPTestCase import DRPTestCase, runTests
from DRP.models import Compound, OrdMolDescriptorValue, NumMolDescriptorValue
from DRP.models import BoolMolDescriptorValue, CatMolDescriptorValue
from DRP.models import OrdMolDescriptor, NumMolDescriptor, BoolMolDescriptor
from DRP.models import CatMolDescriptor, CategoricalDescriptorPermittedValue
from decorators import createsCompound, joinsLabGroup, createsUser
from decorators import createsChemicalClass
from django.core.exceptions import ValidationError
loadTests = unittest.TestLoader().loadTestsFromTestCase


@createsUser('Aslan', 'old_magic')
@joinsLabGroup('Aslan', 'Narnia')
@createsChemicalClass('org', 'Organic')
@createsCompound('EtOH', 682, 'Org', 'Narnia')
class DescriptorsCalced(DRPTestCase):

    """Checks that when a compound is created the descriptors are created."""

    def runTest(self):
        """Run the test."""
        self.assertEqual(1, OrdMolDescriptorValue.objects.count())
        self.assertEqual(1, NumMolDescriptorValue.objects.count())
        self.assertEqual(1, CatMolDescriptorValue.objects.count())
        self.assertEqual(1, BoolMolDescriptorValue.objects.count())


class DoublePluginImport(DRPTestCase):

    """Checks no duplicate key errors on module reload."""

    def runTest(self):
        """Run the test."""
        import DRP.plugins.moldescriptors.example
        reload(DRP.plugins.moldescriptors.example)


class MaxMinValidation(DRPTestCase):

    """Ensures that descriptor minimums and maximums are enforced."""

    def test_ordinal(self):
        """Test ordinal descriptors."""
        with self.assertRaises(ValidationError):
            desc = OrdMolDescriptor(
                heading='heading',
                name='test descriptor',
                calculatorSoftware='test suite',
                calculatorSoftwareVersion=0,
                maximum=3,
                minimum=5
            )
            desc.save()

    def test_ordinal_ok(self):
        """Test a working ordinal descriptor."""
        desc = OrdMolDescriptor(
            heading='heading',
            name='test descriptor',
            calculatorSoftware='test suite',
            calculatorSoftwareVersion=0,
            maximum=5,
            minimum=3)
        desc.save()
        desc.delete()

    def test_numeric(self):
        """Test ordinal descriptors."""
        with self.assertRaises(ValidationError):
            desc = NumMolDescriptor(
                heading='heading',
                name='test descriptor',
                calculatorSoftware='test suite',
                calculatorSoftwareVersion=0,
                maximum=3,
                minimum=5)
            desc.save()

    def test_numeric_ok(self):
        """Test a working numeric descriptor."""
        desc = NumMolDescriptor(
            heading='heading',
            name='test descriptor',
            calculatorSoftware='test suite',
            calculatorSoftwareVersion=0,
            maximum=5,
            minimum=3)
        desc.save()
        desc.delete()

    def test_numeric_max_null(self):
        """Test creating a numeric descriptor with a null max value works."""
        desc = NumMolDescriptor(
            heading='heading',
            name='test descriptor',
            calculatorSoftware='test suite',
            calculatorSoftwareVersion=0,
            maximum=None,
            minimum=3)
        desc.save()
        desc.delete()

    def test_numeric_min_null(self):
        """Test creating a numeric descriptor with a null min value works."""
        desc = NumMolDescriptor(
            heading='heading',
            name='test descriptor',
            calculatorSoftware='test suite',
            calculatorSoftwareVersion=0,
            maximum=5,
            minimum=None)
        desc.save()
        desc.delete()

    def test_numeric_lims_null(self):
        """Test creating a numeric descriptor w/a null min, max value works."""
        desc = NumMolDescriptor(
            heading='heading',
            name='test descriptor',
            calculatorSoftware='test suite',
            calculatorSoftwareVersion=0,
            maximum=None,
            minimum=None)
        desc.save()
        desc.delete()


@createsUser('Aslan', 'old_magic')
@joinsLabGroup('Aslan', 'Narnia')
@createsChemicalClass('Org', 'Organic')
@createsCompound('EtOH', 682, 'Org', 'Narnia')
class MaxMinValueValidation(DRPTestCase):

    """Ensures that descriptor values cannot exceed their prescribed range."""

    def test_numeric_ok(self):
        """Test a working numeric descriptor."""
        desc = NumMolDescriptor(
            heading='heading',
            name='test descriptor',
            calculatorSoftware='test suite',
            calculatorSoftwareVersion=0,
            maximum=5,
            minimum=3)
        desc.save()
        descVal = NumMolDescriptorValue(
            compound=Compound.objects.get(abbrev='EtOH'),
            descriptor=desc,
            value=3)
        descVal.save()
        descVal.delete()
        desc.delete()

    def test_numeric_toohigh(self):
        """Test a descriptor value that is above the maximum."""
        desc = NumMolDescriptor(
            heading='heading',
            name='test descriptor',
            calculatorSoftware='test suite',
            calculatorSoftwareVersion=0,
            maximum=5,
            minimum=3)
        desc.save()
        with self.assertRaises(ValidationError):
            descVal = NumMolDescriptorValue(
                compound=Compound.objects.get(abbrev='EtOH'),
                descriptor=desc,
                value=6)
            descVal.save()
        desc.delete()

    def test_numeric_loolow(self):
        """Test a descriptor value that is below the minimum."""
        desc = NumMolDescriptor(
            heading='heading',
            name='test descriptor',
            calculatorSoftware='test suite',
            calculatorSoftwareVersion=0,
            maximum=5,
            minimum=3)
        desc.save()
        with self.assertRaises(ValidationError):
            descVal = NumMolDescriptorValue(
                compound=Compound.objects.get(abbrev='EtOH'),
                descriptor=desc,
                value=2)
            descVal.save()
        desc.delete()

    def test_numeric_max_null(self):
        """Test creating a numeric descriptor with a null max value works."""
        desc = NumMolDescriptor(
            heading='heading',
            name='test descriptor',
            calculatorSoftware='test suite',
            calculatorSoftwareVersion=0,
            maximum=None,
            minimum=3)
        desc.save()
        descVal = NumMolDescriptorValue(
            compound=Compound.objects.get(abbrev='EtOH'),
            descriptor=desc,
            value=3)
        descVal.save()
        desc.delete()

    def test_numeric_max_null_toolow(self):
        """Test creating a numeric descriptor with a null max value works."""
        desc = NumMolDescriptor(
            heading='heading',
            name='test descriptor',
            calculatorSoftware='test suite',
            calculatorSoftwareVersion=0,
            maximum=None,
            minimum=3)
        desc.save()
        with self.assertRaises(ValidationError):
            descVal = NumMolDescriptorValue(
                compound=Compound.objects.get(abbrev='EtOH'),
                descriptor=desc,
                value=2)
            descVal.save()
        desc.delete()

    def test_numeric_max_null_toohigh(self):
        """Test creating a numeric descriptor with a null max value works."""
        desc = NumMolDescriptor(
            heading='heading',
            name='test descriptor',
            calculatorSoftware='test suite',
            calculatorSoftwareVersion=0,
            maximum=5,
            minimum=None)
        desc.save()
        with self.assertRaises(ValidationError):
            descVal = NumMolDescriptorValue(
                compound=Compound.objects.get(abbrev='EtOH'),
                descriptor=desc,
                value=6)
            descVal.save()
        desc.delete()

    def test_numeric_min_null(self):
        """Test creating a numeric descriptor with a null min value works."""
        desc = NumMolDescriptor(
            heading='heading',
            name='test descriptor',
            calculatorSoftware='test suite',
            calculatorSoftwareVersion=0,
            maximum=5,
            minimum=None)
        desc.save()
        descVal = NumMolDescriptorValue(
            compound=Compound.objects.get(abbrev='EtOH'),
            descriptor=desc,
            value=3)
        descVal.save()
        desc.delete()

    def test_numeric_lims_null(self):
        """Test creating a numeric descriptor w/a null min, max value works."""
        desc = NumMolDescriptor(
            heading='heading',
            name='test descriptor',
            calculatorSoftware='test suite',
            calculatorSoftwareVersion=0,
            maximum=None,
            minimum=None)
        desc.save()
        descVal = NumMolDescriptorValue(
            compound=Compound.objects.get(abbrev='EtOH'),
            descriptor=desc,
            value=3)
        descVal.save()
        descVal.delete()
        desc.delete()

    def test_ordinal_ok(self):
        """Test a working ordinal descriptor."""
        desc = OrdMolDescriptor(
            heading='heading',
            name='test descriptor',
            calculatorSoftware='test suite',
            calculatorSoftwareVersion=0,
            maximum=5,
            minimum=3)
        desc.save()
        descVal = OrdMolDescriptorValue(
            compound=Compound.objects.get(abbrev='EtOH'),
            descriptor=desc,
            value=3)
        descVal.save()
        descVal.delete()
        desc.delete()

    def test_ordinal_toohigh(self):
        """Test a descriptor value that is above the maximum."""
        desc = OrdMolDescriptor(
            heading='heading',
            name='test descriptor',
            calculatorSoftware='test suite',
            calculatorSoftwareVersion=0,
            maximum=5,
            minimum=3)
        desc.save()
        with self.assertRaises(ValidationError):
            descVal = OrdMolDescriptorValue(
                compound=Compound.objects.get(abbrev='EtOH'),
                descriptor=desc,
                value=6)
            descVal.save()
        desc.delete()

    def test_ordinal_loolow(self):
        """Test a descriptor value that is below the minimum."""
        desc = OrdMolDescriptor(
            heading='heading',
            name='test descriptor',
            calculatorSoftware='test suite',
            calculatorSoftwareVersion=0,
            maximum=5,
            minimum=3)
        desc.save()
        with self.assertRaises(ValidationError):
            descVal = OrdMolDescriptorValue(
                compound=Compound.objects.get(abbrev='EtOH'),
                descriptor=desc,
                value=2)
            descVal.save()
        desc.delete()

@createsUser('Aslan', 'old_magic')
@joinsLabGroup('Aslan', 'Narnia')
@createsChemicalClass('Org', 'Organic')
@createsCompound('EtOH', 682, 'Org', 'Narnia')
class CategoricalValidation(DRPTestCase):

    """Ensures that validation for categorical descriptors works."""

    def setUp(self):
        """Set up the test."""
        self.desc = CatMolDescriptor(
            heading='heading',
            name='example',
            calculatorSoftware='test suite',
            calculatorSoftwareVersion=0)
        self.desc.save()
        self.desc2 = CatMolDescriptor(
            heading='heading2',
            name='example2',
            calculatorSoftware='test suite',
            calculatorSoftwareVersion=0)
        self.desc2.save()
        cmdv = CategoricalDescriptorPermittedValue
        self.descPerm = cmdv.objects.get_or_create(
            descriptor=self.desc,
            value='fun')[0]
        self.descPerm.save()
        self.descPerm2 = cmdv.objects.get_or_create(
            descriptor=self.desc,
            value='dull')[0]
        self.descPerm2.save()
        self.descPerm3 = cmdv.objects.get_or_create(
            descriptor=self.desc2,
            value='banana')[0]
        self.descPerm3.save()

    def test_fine(self):
        """Test cases that should work."""
        cmdv = CatMolDescriptorValue(
            descriptor=self.desc,
            compound=Compound.objects.get(CSID=682),
            value=self.descPerm)
        cmdv.save()
        cmdv.delete()

    def test_broken(self):
        """Test broken cases."""
        with self.assertRaises(ValidationError):
            cmdv = CatMolDescriptorValue(
                descriptor=self.desc,
                compound=Compound.objects.get(CSID=682),
                value=self.descPerm3)
            cmdv.save()
        cmdv.delete()

    def tearDown(self):
        """Clean up the test."""
        CatMolDescriptorValue.objects.all().delete()
        self.descPerm.delete()
        self.descPerm2.delete()
        self.descPerm3.delete()
        self.desc.delete()
        self.desc2.delete()

suite = unittest.TestSuite([
    loadTests(DescriptorsCalced),
    loadTests(DoublePluginImport),
    loadTests(MaxMinValidation),
    loadTests(MaxMinValueValidation),
    loadTests(CategoricalValidation)
])

if __name__ == '__main__':
    runTests(suite)
