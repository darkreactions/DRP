#!/usr/bin/env python
'''This module contains tests for teh confirmation page'''

from django.conf import settings
from HttpTest import GetHttpTest, PostHttpTest, GetHttpSessionTest, PostHttpSessionTest
from HttpTest import  OneRedirectionMixin, logsInAs, usesCsrf
from HttpTest import  choosesLabGroup 
from DRP.tests.decorators import joinsLabGroup, createsChemicalClass, signsExampleLicense
from DRP.tests.decorators import createsUser, createsCompound
from DRP.tests import runTests
from django.contrib.auth.models import User
from DRP.models import ConfirmationCode, LabGroup, ChemicalClass, License, Compound
from django.core.urlresolvers import reverse
from uuid import uuid4
import requests
import unittest
from datetime import date, timedelta
from os import path

loadTests = unittest.TestLoader().loadTestsFromTestCase

newCompoundUrl = GetHttpTest.baseUrl + reverse('newCompound')
compoundListUrl = GetHttpTest.baseUrl + reverse('compoundguide', args=['/'])

@logsInAs('Aslan', 'old_magic')
class LicenseRedirect(GetHttpSessionTest, OneRedirectionMixin):
  '''Tests that the request is redirected if a user tries to view the compound add page in without having
  signed an EULA.'''

  url=newCompoundUrl
  testCodes = ['c9e46ba1-cd2a-4080-88b5-97415fa7c484']

  def setUp(self):
    self.license = License(text='This is an example license used in a test', effectiveDate=date.today() - timedelta(1))
    self.license.save()
    super(LicenseRedirect, self).setUp()

  def tearDown(self):
    self.license.delete()

@logsInAs('Aslan', 'old_magic')
@signsExampleLicense('Aslan')
class Lab403Test(GetHttpSessionTest):
  '''Tests that the view returns the special 403 page when trying to look at a
  compound guide without being in a research group
  '''

  url=newCompoundUrl
  status = 403
  testCodes=['91b3d85b-f975-45f1-b0b5-455d475cfa30']

@logsInAs('Aslan', 'old_magic')
@signsExampleLicense('Aslan')
@joinsLabGroup('Aslan', 'narnia')
class CreateCompoundGetTest(GetHttpSessionTest):
  '''Tests that when signed in with full credentials, the
  create view displays'''

  url=newCompoundUrl
  testCodes=['575b31b0-60d1-41d3-86a1-83a8a8b3a7a6', 'd41e5f12-88fd-4494-90fd-96aa84e5beea'] #first one tests for textbox CSID input, second tests correct template

@logsInAs('Aslan', 'old_magic')
@signsExampleLicense('Aslan')
@joinsLabGroup('Aslan', 'Narnia')
@createsChemicalClass('Org', 'Organic')
@usesCsrf
class CreateCompoundRedirTest(PostHttpSessionTest, OneRedirectionMixin):
  '''Tests that the create compound redirection works, and by proxy that the list displays when compounds are present'''

  url=newCompoundUrl
  testCodes=['bf3a3711-b21d-4710-a989-6d1ebc1c9ee9', '7f25b7df-2176-455b-9a68-620af1d52e46']#the first of these tests for correct template, the second tests that the compound table gets displayed 
  _payload = {'abbrev':'etoh', 'name':'ethanol', 'CAS_ID':'64-17-5', 'CSID':'682'}
  
  def setUp(self):
    self.payload['labGroup']=LabGroup.objects.get(title='Narnia').id
    self.payload['chemicalClasses']=[ChemicalClass.objects.get(label='Org').id]
    super(CreateCompoundRedirTest, self).setUp()

@logsInAs('Aslan', 'old_magic')
@signsExampleLicense('Aslan')
@joinsLabGroup('Aslan', 'Narnia')
@usesCsrf
class CreateCompoundRadioTest(PostHttpSessionTest):
  '''Tests for the display of the radio buttons section when presented only with a CSID'''

  url=newCompoundUrl
  testCodes=['1bf53b3a-ddf0-407b-b565-b732e4fa5ddb']#tests for presence of CSID radiobuttons
  _payload={'name':'ethanol'}

@logsInAs('Aslan', 'old_magic')
@signsExampleLicense('Aslan')
@joinsLabGroup('Aslan', 'Narnia')
class NoCompounds(GetHttpSessionTest):
  '''Tests that the empy message is displayed when a group has no compounds'''

  url = compoundListUrl
  testCodes = ['1bf53b3a-ddf0-407b-b565-b732e4fa5ddb']#tests for empty list message

@logsInAs('Aslan', 'old_magic')
@signsExampleLicense('Aslan')
@joinsLabGroup('Aslan', 'Narnia')
@joinsLabGroup('Aslan', 'Stone Table')
class ManyGroupsRedirect(GetHttpSessionTest, OneRedirectionMixin):
  '''Tests that a user with many lab groups but no session data for a lab group gets redirected. Tests the display of the lab group selection template by proxy.'''

  url = compoundListUrl
  testCodes = ['82ab2a5b-d337-4579-89d4-621cf2ce07ea']

@logsInAs('Aslan', 'old_magic')
@signsExampleLicense('Aslan')
@joinsLabGroup('Aslan', 'Narnia')
@joinsLabGroup('Aslan', 'Stone Table')
@choosesLabGroup('Aslan', 'Narnia')
class ManyLabGroupsDisplays(GetHttpSessionTest): 
  '''Tests that a user with many lab groups with session data for a lab group does not get redirected'''

  url=compoundListUrl
  testCodes = ['bf3a3711-b21d-4710-a989-6d1ebc1c9ee9']

@logsInAs('Aslan', 'old_magic')
@signsExampleLicense('Aslan')
@joinsLabGroup('Aslan', 'Narnia')
@joinsLabGroup('Aslan', 'Stone table')
@usesCsrf
class LabGroupSelectionRedirect(PostHttpSessionTest, OneRedirectionMixin):
  '''tests for the redirection after the choice of lab group has been made'''

  url = PostHttpSessionTest.baseUrl + reverse('selectGroup')
  testCodes = ['bf3a3711-b21d-4710-a989-6d1ebc1c9ee9']
  _params={'next':compoundListUrl}

  def setUp(self, *args, **kwargs):
    self.payload['labGroup'] = LabGroup.objects.get(title='Narnia').id
    super(LabGroupSelectionRedirect, self).setUp(*args, **kwargs)

@logsInAs('Aslan', 'old_magic')
@signsExampleLicense('Aslan')
@joinsLabGroup('Aslan', 'Narnia')
@createsChemicalClass('Org', 'Organic')
@createsCompound('EtOH', 682, 'Org', 'Narnia')
class GetCompoundForEditing(GetHttpSessionTest):
  '''Tests that fetching a compound for editing works'''

  testCodes = ['7d3763bc-c7d0-4102-a036-8c184263fe21']

  def setUp(self):
    self.url = self.baseUrl + reverse('editCompound', args=[Compound.objects.get(abbrev='EtOH').pk])
    super(GetCompoundForEditing, self).setUp()

@logsInAs('Aslan', 'old_magic')
@createsUser('White Witch', 'new_magic')
@signsExampleLicense('Aslan')
@joinsLabGroup('Aslan', 'Narnia')
@joinsLabGroup('White Witch', 'stone table')
@createsChemicalClass('Org', 'Organic')
@createsCompound('EtOH', 682, 'Org', 'Narnia')
@createsCompound('Pyr', 8904, 'Org', 'stone table')
class GetNotMyCompoundForEditing(GetHttpSessionTest):
  '''Tests that fetching someone elses compound returns a 404''' 

  status=404

  def setUp(self):
    self.url = self.baseUrl + reverse('editCompound', args=[Compound.objects.get(abbrev='Pyr').pk])
    super(GetNotMyCompoundForEditing, self).setUp()

@logsInAs('Aslan', 'old_magic')
@signsExampleLicense('Aslan')
@joinsLabGroup('Aslan', 'Narnia')
@createsChemicalClass('Org', 'Organic')
@createsCompound('EtOH', 682, 'Ethanol', 'Narnia', custom=True)
class GetCustomCompound403(GetHttpSessionTest):
  '''Tests that fetching a compound with the custom flag for editing returns a 403'''

  status=403 

  def setUp(self):
    self.url = self.baseUrl + reverse('editCompound', args=[Compound.objects.get(abbrev='EtOH').pk])
    super(GetCustomCompound403, self).setUp()

@logsInAs('Aslan', 'old_magic')
@signsExampleLicense('Aslan')
@joinsLabGroup('Aslan', 'Narnia')
class GetCompoundUpload(GetHttpSessionTest):
  '''Tests that GETing the compound upload page works'''

  url = GetHttpSessionTest.baseUrl + reverse('uploadcompoundcsv')
  testCodes = ["d16ff9ed-752c-405e-bc4b-cae3a27dd7b2"] 

@logsInAs('Aslan', 'old_magic')
@signsExampleLicense('Aslan')
@joinsLabGroup('Aslan', 'Narnia')
@usesCsrf
class PostCompoundUpload(PostHttpSessionTest, OneRedirectionMixin):

  url = PostHttpSessionTest.baseUrl + reverse('uploadcompoundcsv')
  testCodes = ["bf3a3711-b21d-4710-a989-6d1ebc1c9ee9"] 
  
  def setUp(self):
    self.payload['labGroup'] = LabGroup.objects.get(title="Narnia").pk
    with open(path.join(settings.APP_DIR, 'tests', 'resource', 'compound_spread_test1.csv'), 'rb') as f:
      self.files['csv']=('compound_test1.csv',f.read(), 'text/csv')
    super(PostCompoundUpload, self).setUp()

# TODO XXX PHIL_FIX_AFTER_CONTEXT_PROCESSOR
suite = unittest.TestSuite([
  loadTests(LicenseRedirect),
  loadTests(Lab403Test),
  loadTests(CreateCompoundGetTest),
  #loadTests(CreateCompoundRedirTest),
  loadTests(CreateCompoundRadioTest),
  loadTests(NoCompounds),
  loadTests(ManyGroupsRedirect),
  loadTests(ManyLabGroupsDisplays),
  loadTests(LabGroupSelectionRedirect),
  loadTests(GetCompoundForEditing),
  loadTests(GetNotMyCompoundForEditing),
  loadTests(GetCustomCompound403),
  loadTests(GetCompoundUpload),
  #loadTests(PostCompoundUpload)
])

if __name__=='__main__':
  runTests(suite)
