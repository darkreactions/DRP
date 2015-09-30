#!/usr/bin/env python
'''The unit test for the compound editing form.
  These tests assume that presence tests for the form fields work as expected
'''

import unittest 
from BaseFormTest import BaseFormTest
from DRP.forms import CompoundEditForm
from DRP.models import LabGroup, ChemicalClass, Compound, LabGroup
from django.conf import settings
from DRP.tests.decorators import createsCompound, createsUser, joinsLabGroup, createsChemicalClass 
from DRP.tests import runTests
loadTests = unittest.TestLoader().loadTestsFromTestCase

@createsUser('Aslan', 'old_magic')
@joinsLabGroup('Aslan', 'Narnia')
@createsChemicalClass('Org', 'Organic Reagent')
@createsCompound('EtOH', 682, 'Org', 'Narnia')
class CorrectSynonym(BaseFormTest):
  '''Tests that the form validates with a correct synonym for the compound submitted'''

  def test_validation(self):
    self.validationSucceeds()

  def setUpFormData(self):
    self.formData = {'name':'Ethanol', 'abbrev':'ban'}
    self.formData['chemicalClasses'] = [c.id for c in ChemicalClass.objects.filter(label='Org')]

  def setUp(self):
    "Instantiates the form"
    super(CorrectSynonym, self).setUp()
    self.form = CompoundEditForm(data=self.formData,instance=Compound.objects.get(abbrev='EtOH', labGroup=LabGroup.objects.get(title="Narnia")))

@createsUser('Aslan', 'old_magic')
@joinsLabGroup('Aslan', 'Narnia')
@createsChemicalClass('Org', 'Organic Reagent')
@createsCompound('EtOH', 682, 'Org', 'Narnia')
class IncorrectSynonym(BaseFormTest):
  '''Tests that the form fails to validate when provided with an incorrect synonym'''
  
  def setUpFormData(self):
    self.formData = {'name':'Pyrazine', 'abbrev':'ban'}
    self.formData['chemicalClasses'] = [c.id for c in ChemicalClass.objects.filter(label='Org')]

  def setUp(self):
    "Instantiates the form"
    super(IncorrectSynonym, self).setUp()
    self.form = CompoundEditForm(data=self.formData, instance=Compound.objects.get(abbrev='EtOH', labGroup=LabGroup.objects.get(title="Narnia")))

suite = unittest.TestSuite([
          loadTests(IncorrectSynonym),
          loadTests(CorrectSynonym)
          ])

if __name__=='__main__':
  runTests(suite)
