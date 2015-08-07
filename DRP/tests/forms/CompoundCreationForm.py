#!/usr/bin/env python
'''The unit test for the LabGroupForm class.
  These tests assume that presence tests for teh form fields work as expected
  This test suite assumes that Django's check_password functions as
  expected.
  This test suite assumes that Django's save methods on forms work for
  trivial cases.
'''

import unittest 
from DRP.tests import DRPTestCase
from DRP.forms import CompoundForm
from DRP.models import LabGroup, ChemicalClass
from django.conf import settings
from django.contrib.auth.models import User
loadTests = unittest.TestLoader().loadTestsFromTestCase

class NoLabExists(DRPTestCase):
  '''Tests that the form doesn't validate when there are no lab groups'''
  

  def setUpFormData(self):
    self.formData = {'labGroup':'5', 'abbrev':'etoh', 'name':'ethanol', 'CAS_ID':'64-17-5', 'CSID':'682'}
    self.formData['chemicalClass'] = [str(self.chemicalClass.id)]

  def setUp(self):
    '''Creates a user and a chemical class, then a form'''
    self.user = User.objects.create_user('Aslan', 'old_magic')
    self.user.save()
    self.chemicalClass = ChemicalClass(label='Solv', description='Common Solvent')
    self.chemicalClass.save()    
    self.setUpFormData()
    self.form = CompoundForm(self.user, self.formData)

  def test_validation(self):
    self.validationFails()

  def validationFails(self):
    '''a test for cases where the form should not validate'''
    self.assertFalse(self.form.is_valid(), 'Form should have failed... submitted data: {0}'.format(self.formData))

  def validationSucceeds(self):
    '''Test that the form does validates'''
    valid = self.form.is_valid()
    errString = ''
    for e, m in self.form.errors.items():
      errString += '{0}: {1}\n'.format(e, m)
    errString += ' Submitted data = {0}'.format(self.formData)
    self.assertTrue(valid, errString)

  def tearDown(self):
    self.chemicalClass.delete()
    self.user.delete()

class NoLabForUser(NoLabExists):
  '''Tests that the form doesn't validate when there are lab groups but the user is not a member'''
  
  def setUpFormData(self):
    super(NoLabForUser, self).setUpFormData()
    self.formData['labGroup']=str(self.labGroup.id)

  def setUp(self):
    self.labGroup = LabGroup.objects.makeLabGroup(title="LegacyPassTest1", address='1, war drobe, Narnia', email='aslan@example.com', access_code='old_magic')
    self.labGroup.save()
    super(NoLabForUser, self).setUp()

  def tearDown(self):
    self.labGroup.delete()
    super(NoLabForUser, self).tearDown()

class LabForUser(NoLabForUser):
  '''Tests that the form validates when the user is a member of the lab. Also checks that the save method works correctly.'''

  compound = None
  
  def setUpFormData(self):
    super(LabForUser, self).setUpFormData()
    self.labGroup.users.add(self.user)
    self.labGroup.save()
  
  def test_validation(self):
    self.validationSucceeds()

  def test_saving(self):
    if self.form.is_valid():
      self.compound = self.form.save()
      self.assertIsNotNone(self.compound.id)

  def tearDown(self):
    super(LabForUser, self).tearDown()
    if self.compound is not None:
      self.compound.delete()

class LabForBadUser(LabForUser):
  '''Tests that the form doesn't validate when the user doesn't choose a lab group'''

  def setUpFormData(self):
    super(LabForBadUser, self).setUpFormData()
    self.formData['labGroup'] = ''

  def test_validation(self):
    self.validationFails()

class NoCAS(LabForUser):
  '''Tests that the form validates with no CAS number provided'''

  def setUpFormData(self):
    super(NoCAS, self).setUpFormData()
    self.formData['CAS_ID'] = ''
  
class InconsistentCASName(LabForUser):
  '''Tests that the form does not validate when the CAS number and name are not consistent'''

  def setUpFormData(self):
    super(InconsistentCASName, self).setUpFormData()
    self.formData['CAS_ID'] = '290-37-9'
    self.formData['CSID'] = '8904' #consistent values for pyrazine.

  def test_validation(self):
    self.validationFails()

class InconsistentNameCSID(LabForUser):
  '''Tests that the form does not validate when the name and CSID are not consistent''' 
  
  def setUpFormData(self):
    super(InconsistentNameCSID, self).setUpFormData()
    self.formData['name'] = 'Pyrazine'
  
  def test_validation(self):
    self.validationFails()

class InconsistentCSIDCAS(LabForUser):
  '''Tests that the form does not validate when the CAS number and CSID are not consistent'''

  def setUpFormData(self):
    super(InconsistentCSIDCAS, self).setUpFormData()
    self.formData['CAS_ID'] = '290-37-9'

  def test_validation(self):
    self.validationFails()

suite = unittest.TestSuite([
          loadTests(NoLabExists),
          loadTests(NoLabForUser),
          loadTests(LabForUser),
          loadTests(LabForBadUser),
          loadTests(NoCAS),
          loadTests(InconsistentCASName),
          loadTests(InconsistentNameCSID),
          loadTests(InconsistentCSIDCAS)
          ])

if __name__=='__main__':
  unittest.TextTestRunner(verbosity=2).run(suite)
