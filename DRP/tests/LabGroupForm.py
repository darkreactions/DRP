#!/usr/bin/env python
'''The unit test for the LabGroupForm class.
  This test suite assumes that Django's inbuilt email testing is already
  conformant to the relevant standards.
  This test suite assumes that Django's check_password functions as
  expected.
  This test suite assumes that Django's save methods on forms work for
  trivial cases.
'''

import unittest 
import TestConfig
from DRP.forms import LabGroupForm
from DRP.models import LabGroup
from DRP.settings import LAB_GROUP_HASH_SALT 
from django.contrib.auth.hashers import check_password, make_password
from django.contrib.auth.models import User
from copy import copy

loadTests = unittest.TestLoader.loadTestsFromTestCase

class LegacyPassword(unittest.TestCase):
  '''Tests for a case where a legacy password is set but no new password has been entered'''

  def setUp(self):
    'instantiate the labgroup'''
    self.labGroup = LabGroup(title="LegacyPassTest1", address='1, war drobe, Narnia', email='aslan@example.com', legacy_access_code='old_magic')
    formData = {'title':'LegacyPasstest1', 'address':'1, war drobe', 'email':'aslan@example.com', 'users':[]}
    self.form = LabGroupForm(formData, instance=self.labGroup)
    self.labGroup.save()

  def test_validation(self):
    valid = self.form.is_valid()
    errString = ''
    for e, m in self.form.errors.items():
      errString += '{0}: {1}\n'.format(e, m)
    self.assertTrue(valid, errString)

  def test_databaseChange(self):
    self.form.save()
    self.labGroup = LabGroup.objects.get(pk=self.labGroup.id)
    self.assertEqual(self.labGroup.legacy_access_code, '')
    self.assertTrue(check_password('old_magic', self.labGroup.access_code))

  def tearDown(self):
    self.labGroup.delete()

class LegacyPassword2(unittest.TestCase):
  '''Tests for a case where a legacy password is set and a new password has been entered'''

  def setUp(self):
    '''instantiate the labgroup'''
    self.labGroup = LabGroup(title="LegacyPassTest2", address='2, war drobe, Narnia', email='whitewitch@example.com', legacy_access_code='old_magic')
    self.labGroup.save()
    formData = {'title':'LegacyPasstest2', 'address':'2, war drobe', 'email':'whitewitch@example.com', 'accessCode':'new_magic'}
    self.form = LabGroupForm(formData, instance=self.labGroup)

  def test_validation(self):
    self.assertTrue(self.form.is_valid())

  def test_databaseChange(self):
    self.form.save()
    self.labGroup= LabGroup.objects.get(pk=self.labGroup.id)
    self.assertEqual(self.labGroup.legacy_access_code, '')
    self.assertTrue(check_password('new_magic', self.labGroup.access_code))

  def tearDown(self):
    self.labGroup.delete()

class CreateNew(unittest.TestCase):

  def setUp(self):
    self.user1 = User(first_name='Aslan', password=make_password('old_magic'), username="Aslan", email="aslan@example.com")
    self.user2 = User(first_name='White', last_name='Witch', password=make_password('new_magic'), username="whitewitch", email="whitewitch@example.com")
    self.user1.save()
    self.user2.save()
    self.users = [self.user1, self.user2]
    self.formData = {'title':'CreationTest', 'address':'war drobe', 'email':'whitewitch@example.com', 'accessCode':'turkishdelight'}

  def test_missingTitle(self):
    formData = copy(self.formData)
    del formData['title']
    form = LabGroupForm(formData)
    self.assertFalse(form.is_valid())

  def test_missingAddress(self):
    formData = copy(self.formData)
    del formData['address']
    form = LabGroupForm(formData)
    self.assertFalse(form.is_valid())

  def test_missingEmail(self):
    formData = copy(self.formData)
    del formData['email']
    form = LabGroupForm(formData)
    self.assertFalse(form.is_valid())

  def test_missingPass(self):
    formData = copy(self.formData)
    del formData['accessCode']
    form = LabGroupForm(formData)
    self.assertFalse(form.is_valid())

  def completeData(self):
    form = LabGroupForm(self.formData)
    self.assertTrue(form.is_valid())

  def completeDataAndUsers(self):
    formData = copy(self.formData)
    formData['users'] = [u.id for u in self.users]
    form = LabGroupForm()
    self.assertTrue(form.is_valid())

  def tearDown(self):
    self.user1.delete()
    self.user2.delete()

class DuplicateUniqueValues(unittest.TestCase): 

  def setUp(self):
    self.labGroup = LabGroup(title="DuplicateValues", address='2, war drobe, Narnia', email='whitewitch@example.com', legacy_access_code='old_magic')
    self.labGroup.save()
  
  def test_duplicateTitle(self):
    form = LabGroupForm({'title':'DuplicateValues', 'address':'2, war drobe, Narnia', 'email':'whitewitch@example.com', 'accessCode':'old_magic'})
    self.assertFalse(form.is_valid())
    
  def tearDown(self):
    self.labGroup.delete()

def suite():
  #This function should be adjusted to contain the loadTests() function enacted on each test case.
  return unittest.TestSuite([
          loadTests(LegacyPassword),
          loadTests(LegacyPassword2),
          loadTests(CreateNew),
          loadTests(DuplicateUniqueValues)
          ])

if __name__ == '__main__':
  #Runs the test- a good way to check that this particular test set works without having to run all the tests.
  unittest.main()
