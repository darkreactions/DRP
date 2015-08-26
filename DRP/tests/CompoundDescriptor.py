#!/usr/bin/env python
'''A module containing tests for molecular descriptor classes'''

import unittest
from DRPTestCase import DRPTestCase, runTests
from DRP.models import Compound, OrdMolDescriptorValue, NumMolDescriptorValue, BoolMolDescriptorValue, CatMolDescriptorValue
from DRP.models import OrdMolDescriptor, NumMolDescriptor, BoolMolDescriptor, CatMolDescriptor, CatMolDescriptorPermitted
from decorators import createsCompound, joinsLabGroup, createsUser, createsChemicalClass
from django.core.exceptions import ValidationError
loadTests = unittest.TestLoader().loadTestsFromTestCase

@createsUser('Aslan', 'old_magic')
@joinsLabGroup('Aslan', 'Narnia')
@createsChemicalClass('org', 'Organic')
@createsCompound('EtOH', 682, 'Org', 'Narnia') 
class DescriptorsCalced(DRPTestCase):
  '''Checks that when a compound is created the right number of descriptors are created'''

  def runTest(self):
    self.assertEqual(1, OrdMolDescriptorValue.objects.count())
    self.assertEqual(1, NumMolDescriptorValue.objects.count())
    self.assertEqual(1, CatMolDescriptorValue.objects.count())
    self.assertEqual(1, BoolMolDescriptorValue.objects.count())
      
class DoublePluginImport(DRPTestCase):
  '''makes sure that when the plugin is reloaded we don't wind up with duplicate key errors'''

  def runTest(self):
    import DRP.plugins.moldescriptors.example
    reload(DRP.plugins.moldescriptors.example)

class MaxMinValidation(DRPTestCase):
  '''Ensures that Descriptors with minimum values higher than the maximums cannot be created'''

  def test_ordinal(self):
    '''tests ordinal descriptors'''
    with self.assertRaises(ValidationError):
      desc = OrdMolDescriptor(heading='heading', name='test descriptor', calculatorSoftware='test suite', calculatorSoftwareVersion=0, maximum=3, minimum=5) 
      desc.save()

  def test_ordinal_ok(self):
    '''tests a working ordinal descriptor'''
    desc = OrdMolDescriptor(heading='heading', name='test descriptor', calculatorSoftware='test suite', calculatorSoftwareVersion=0, maximum=5, minimum=3) 
    desc.save()
    desc.delete()
  
  def test_ordinal_max_null(self):
    '''tests that creating a ordinal descriptor with a null max value works'''
    desc = OrdMolDescriptor(heading='heading', name='test descriptor', calculatorSoftware='test suite', calculatorSoftwareVersion=0, maximum=None, minimum=3) 
    desc.save()
    desc.delete()

  def test_ordinal_min_null(self):
    '''tests that creating a ordinal descriptor with a null min value works'''
    desc = OrdMolDescriptor(heading='heading', name='test descriptor', calculatorSoftware='test suite', calculatorSoftwareVersion=0, maximum=5, minimum=None) 
    desc.save()
    desc.delete()

  def test_ordinal_lims_null(self):
    desc = OrdMolDescriptor(heading='heading', name='test descriptor', calculatorSoftware='test suite', calculatorSoftwareVersion=0, maximum=None, minimum=None) 
    desc.save()
    desc.delete()
    '''tests that creating a ordinal descriptor with a null min and max value works'''
  
  def test_numeric(self): 
    '''tests ordinal descriptors'''
    with self.assertRaises(ValidationError):
      desc = NumMolDescriptor(heading='heading', name='test descriptor', calculatorSoftware='test suite', calculatorSoftwareVersion=0, maximum=3, minimum=5) 
      desc.save()

  def test_numeric_ok(self): 
    '''tests a working numeric descriptor'''
    desc = NumMolDescriptor(heading='heading', name='test descriptor', calculatorSoftware='test suite', calculatorSoftwareVersion=0, maximum=5, minimum=3) 
    desc.save()
    desc.delete()

  def test_numeric_max_null(self):
    '''tests that creating a numeric descriptor with a null max value works'''
    desc = NumMolDescriptor(heading='heading', name='test descriptor', calculatorSoftware='test suite', calculatorSoftwareVersion=0, maximum=None, minimum=3) 
    desc.save()
    desc.delete()

  def test_numeric_min_null(self):
    '''tests that creating a numeric descriptor with a null min value works'''
    desc = NumMolDescriptor(heading='heading', name='test descriptor', calculatorSoftware='test suite', calculatorSoftwareVersion=0, maximum=5, minimum=None) 
    desc.save()
    desc.delete()

  def test_numeric_lims_null(self):
    '''tests that creating a numeric descriptor with a null min and max value works'''
    desc = NumMolDescriptor(heading='heading', name='test descriptor', calculatorSoftware='test suite', calculatorSoftwareVersion=0, maximum=None, minimum=None) 
    desc.save()
    desc.delete()

@createsUser('Aslan', 'old_magic')
@joinsLabGroup('Aslan', 'Narnia')
@createsChemicalClass('Org', 'Organic')
@createsCompound('EtOH', 682, 'Org', 'Narnia') 
class MaxMinValueValidation(DRPTestCase):
  '''Ensures that descriptor values cannot exceed their prescribed range, and that null maxs and mins don't problematicise this.'''

  def test_numeric_ok(self): 
    '''tests a working numeric descriptor'''
    desc = NumMolDescriptor(heading='heading', name='test descriptor', calculatorSoftware='test suite', calculatorSoftwareVersion=0, maximum=5, minimum=3) 
    desc.save()
    descVal = NumMolDescriptorValue(compound=Compound.objects.get(abbrev='EtOH'), descriptor=desc, value=3)
    descVal.save()
    descVal.delete()
    desc.delete()

  def test_numeric_toohigh(self):
    '''tests a descriptor value that is above the maximum'''
    desc = NumMolDescriptor(heading='heading', name='test descriptor', calculatorSoftware='test suite', calculatorSoftwareVersion=0, maximum=5, minimum=3) 
    desc.save()
    with self.assertRaises(ValidationError):
      descVal = NumMolDescriptorValue(compound=Compound.objects.get(abbrev='EtOH'), descriptor=desc, value=6)
      descVal.save()
    desc.delete()

  def test_numeric_loolow(self):
    desc = NumMolDescriptor(heading='heading', name='test descriptor', calculatorSoftware='test suite', calculatorSoftwareVersion=0, maximum=5, minimum=3) 
    desc.save()
    with self.assertRaises(ValidationError):
      descVal = NumMolDescriptorValue(compound=Compound.objects.get(abbrev='EtOH'), descriptor=desc, value=2)
      descVal.save()
    desc.delete()

  def test_numeric_max_null(self):
    '''tests that creating a numeric descriptor with a null max value works'''
    desc = NumMolDescriptor(heading='heading', name='test descriptor', calculatorSoftware='test suite', calculatorSoftwareVersion=0, maximum=None, minimum=3) 
    desc.save()
    descVal = NumMolDescriptorValue(compound=Compound.objects.get(abbrev='EtOH'), descriptor=desc, value=3)
    descVal.save()
    desc.delete()

  def test_numeric_max_null_toolow(self):
    '''tests that creating a numeric descriptor with a null max value works'''
    desc = NumMolDescriptor(heading='heading', name='test descriptor', calculatorSoftware='test suite', calculatorSoftwareVersion=0, maximum=None, minimum=3) 
    desc.save()
    with self.assertRaises(ValidationError):
      descVal = NumMolDescriptorValue(compound=Compound.objects.get(abbrev='EtOH'), descriptor=desc, value=2)
      descVal.save()
    desc.delete()
  
  def test_numeric_max_null_toohigh(self):
    '''tests that creating a numeric descriptor with a null max value works'''
    desc = NumMolDescriptor(heading='heading', name='test descriptor', calculatorSoftware='test suite', calculatorSoftwareVersion=0, maximum=5, minimum=None) 
    desc.save()
    with self.assertRaises(ValidationError):
      descVal = NumMolDescriptorValue(compound=Compound.objects.get(abbrev='EtOH'), descriptor=desc, value=6)
      descVal.save()
    desc.delete()

  def test_numeric_min_null(self):
    '''tests that creating a numeric descriptor with a null min value works'''
    desc = NumMolDescriptor(heading='heading', name='test descriptor', calculatorSoftware='test suite', calculatorSoftwareVersion=0, maximum=5, minimum=None) 
    desc.save()
    descVal = NumMolDescriptorValue(compound=Compound.objects.get(abbrev='EtOH'), descriptor=desc, value=3)
    descVal.save()
    desc.delete()

  def test_numeric_lims_null(self):
    '''tests that creating a numeric descriptor with a null min and max value works'''
    desc = NumMolDescriptor(heading='heading', name='test descriptor', calculatorSoftware='test suite', calculatorSoftwareVersion=0, maximum=None, minimum=None) 
    desc.save()
    descVal = NumMolDescriptorValue(compound=Compound.objects.get(abbrev='EtOH'), descriptor=desc, value=3)
    descVal.save()
    descVal.delete()
    desc.delete()

  def test_ordinal_ok(self): 
    '''tests a working ordinal descriptor'''
    desc = OrdMolDescriptor(heading='heading', name='test descriptor', calculatorSoftware='test suite', calculatorSoftwareVersion=0, maximum=5, minimum=3) 
    desc.save()
    descVal = OrdMolDescriptorValue(compound=Compound.objects.get(abbrev='EtOH'), descriptor=desc, value=3)
    descVal.save()
    descVal.delete()
    desc.delete()

  def test_ordinal_toohigh(self):
    '''tests a descriptor value that is above the maximum'''
    desc = OrdMolDescriptor(heading='heading', name='test descriptor', calculatorSoftware='test suite', calculatorSoftwareVersion=0, maximum=5, minimum=3) 
    desc.save()
    with self.assertRaises(ValidationError):
      descVal = OrdMolDescriptorValue(compound=Compound.objects.get(abbrev='EtOH'), descriptor=desc, value=6)
      descVal.save()
    desc.delete()

  def test_ordinal_loolow(self):
    desc = OrdMolDescriptor(heading='heading', name='test descriptor', calculatorSoftware='test suite', calculatorSoftwareVersion=0, maximum=5, minimum=3) 
    desc.save()
    with self.assertRaises(ValidationError):
      descVal = OrdMolDescriptorValue(compound=Compound.objects.get(abbrev='EtOH'), descriptor=desc, value=2)
      descVal.save()
    desc.delete()

  def test_ordinal_max_null(self):
    '''tests that creating a ordinal descriptor with a null max value works'''
    desc = OrdMolDescriptor(heading='heading', name='test descriptor', calculatorSoftware='test suite', calculatorSoftwareVersion=0, maximum=None, minimum=3) 
    desc.save()
    descVal = OrdMolDescriptorValue(compound=Compound.objects.get(abbrev='EtOH'), descriptor=desc, value=3)
    descVal.save()
    desc.delete()

  def test_ordinal_max_null_toolow(self):
    '''tests that creating a ordinal descriptor with a null max value works'''
    desc = OrdMolDescriptor(heading='heading', name='test descriptor', calculatorSoftware='test suite', calculatorSoftwareVersion=0, maximum=None, minimum=3) 
    desc.save()
    with self.assertRaises(ValidationError):
      descVal = OrdMolDescriptorValue(compound=Compound.objects.get(abbrev='EtOH'), descriptor=desc, value=2)
      descVal.save()
    desc.delete()
  
  def test_ordinal_max_null_toohigh(self):
    '''tests that creating a ordinal descriptor with a null max value works'''
    desc = OrdMolDescriptor(heading='heading', name='test descriptor', calculatorSoftware='test suite', calculatorSoftwareVersion=0, maximum=5, minimum=None) 
    desc.save()
    with self.assertRaises(ValidationError):
      descVal = OrdMolDescriptorValue(compound=Compound.objects.get(abbrev='EtOH'), descriptor=desc, value=6)
      descVal.save()
    desc.delete()

  def test_ordinal_min_null(self):
    '''tests that creating a ordinal descriptor with a null min value works'''
    desc = OrdMolDescriptor(heading='heading', name='test descriptor', calculatorSoftware='test suite', calculatorSoftwareVersion=0, maximum=5, minimum=None) 
    desc.save()
    descVal = OrdMolDescriptorValue(compound=Compound.objects.get(abbrev='EtOH'), descriptor=desc, value=3)
    descVal.save()
    desc.delete()

  def test_ordinal_lims_null(self):
    '''tests that creating a ordinal descriptor with a null min and max value works'''
    desc = OrdMolDescriptor(heading='heading', name='test descriptor', calculatorSoftware='test suite', calculatorSoftwareVersion=0, maximum=None, minimum=None) 
    desc.save()
    descVal = OrdMolDescriptorValue(compound=Compound.objects.get(abbrev='EtOH'), descriptor=desc, value=3)
    descVal.save()
    desc.delete()

@createsUser('Aslan', 'old_magic')
@joinsLabGroup('Aslan', 'Narnia')
@createsChemicalClass('Org', 'Organic')
@createsCompound('EtOH', 682, 'Org', 'Narnia') 
class CategoricalValidation(DRPTestCase):
  '''Ensures that validation for categorical descriptors works'''

  def setUp(self):
    self.desc = CatMolDescriptor(heading='heading', name='example', calculatorSoftware='test suite', calculatorSoftwareVersion=0)
    self.desc.save()
    self.desc2 = CatMolDescriptor(heading='heading2', name='example2', calculatorSoftware='test suite', calculatorSoftwareVersion=0)
    self.desc2.save()
    self.descPerm = CatMolDescriptorPermitted.objects.get_or_create(descriptor=self.desc, value='fun')[0]
    self.descPerm.save()
    self.descPerm2 = CatMolDescriptorPermitted.objects.get_or_create(descriptor=self.desc, value='dull')[0]
    self.descPerm2.save()
    self.descPerm3 = CatMolDescriptorPermitted.objects.get_or_create(descriptor=self.desc2, value='banana')[0]
    self.descPerm3.save()
    
  def test_fine(self):
    cmdv = CatMolDescriptorValue(descriptor=self.desc, compound=Compound.objects.get(CSID=682), value=self.descPerm)
    cmdv.save()
    cmdv.delete()
     
  def test_broken(self):
    with self.assertRaises(ValidationError):
      cmdv = CatMolDescriptorValue(descriptor=self.desc, compound=Compound.objects.get(CSID=682), value=self.descPerm3)
      cmdv.save()

  def tearDown(self):
    self.desc.delete()
    self.desc2.delete()

suite = unittest.TestSuite([
          loadTests(DescriptorsCalced),
          loadTests(DoublePluginImport),
          loadTests(MaxMinValidation),
          loadTests(MaxMinValueValidation),
          loadTests(CategoricalValidation)
          ])

if __name__=='__main__':
  runTests(suite)
