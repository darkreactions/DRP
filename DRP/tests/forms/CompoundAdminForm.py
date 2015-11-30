#!/usr/bin/env python
'''Test classes for the compound administration form'''

import unittest 
from DRP.forms import CompoundAdminForm
from BaseFormTest import BaseFormTest
from DRP.models import LabGroup, ChemicalClass
from django.conf import settings
from django.contrib.auth.models import User

loadTests = unittest.TestLoader().loadTestsFromTestCase

class CreateTest(BaseFormTest):
  '''Tests that the form saves data correctly'''

  def setUpFormData(self):
    self.formData = {'labGroup':self.labGroup.id, 'chemicalClasses':[self.chemicalClass.id], 'abbrev':'etoh', 'name':'ethanol', 'CAS_ID':'64-17-5',
                    'CSID':'682', 'INCHI':r'1S/C2H6O/c1-2-3/h3H,2H2,1H3', 'smiles':'CCO'}

  def setUp(self):
    '''Creates a user and a chemical class, then a form'''
    self.user = User.objects.create_user('Aslan', 'old_magic')
    self.user.save()
    self.chemicalClass = ChemicalClass(label='Solv', description='Common Solvent')
    self.chemicalClass.save()    
    self.labGroup = LabGroup.objects.makeLabGroup(title="LegacyPassTest1", address='1, war drobe, Narnia', email='aslan@example.com', access_code='old_magic')
    self.labGroup.save()
    self.setUpFormData()
    self.form = CompoundAdminForm(self.formData)

  def test_validation(self):
    self.validationSucceeds()
  
  def test_saving(self):
    if self.form.is_valid():
      compound = self.form.save()
      self.assertIsNotNone(compound.id)
      self.assertEqual(compound.custom, True)
      compound.delete()

  def tearDown(self):
    self.user.delete()
    self.chemicalClass.delete()
    self.labGroup.delete()

suite = unittest.TestSuite([
  loadTests(CreateTest)
])

if __name__=='__main__':
  unittest.TextTestRunner(verbosity=2).run(suite)
