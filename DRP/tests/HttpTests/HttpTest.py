'''A module containing base classes and decorators for HttpTests'''

import requests
from django.conf import settings
from django.core.urlresolvers import reverse
from django.contrib.auth.models import User
from abc import ABCMeta
from DRP.tests import DRPTestCase
from DRP.models import License, LicenseAgreement, LabGroup, ChemicalClass
import json
from datetime import date, timedelta

class GetHttpTest(DRPTestCase):

  baseUrl = 'http://' + settings.SERVER_NAME
  url = baseUrl
  testCodes = []
  _params = {}
  '''any GET params to be added to the reuqest.'''
  status = 200
  '''The expected status code for this test case'''

  def __init__(self, *args, **kwargs):
    super(GetHttpTest, self).__init__(*args, **kwargs)
    self.params = self._params.copy()

  def setUp(self):
    '''Sets up the test by requesting the home page uri'''
    self.response = requests.get(self.url, params=self.params)

  @staticmethod
  def constructFailureMessage(message):
    '''Constructs failure messages by parsing the output from the w3c validator'''
    m = message['type'] + ':'
    if 'subtype' in message.keys():
      m += message['subtype'] + ':'
    if 'message' in message.keys():
      m += message['message']
    if 'line' in message.keys():
      m += 'on line {0}'.format(message['line'])
    if 'column' in message.keys():
      m += ', column {0}'.format(message['column'])
    return m +'\n'

  def validate(self):
    self.validationResponse = requests.post(settings.EXTERNAL_HTML_VALIDATOR, headers={'content-type':'text/html; charset=utf-8'}, data=self.response.text, params={'out':'json'})

  def test_Status(self):
    '''Checks that the http response code is the expected value'''
    self.assertEqual(self.response.status_code, self.status, 'Url {0} returns code {1}. Page content follows:\n\n{2}'.format(self.url, self.response.status_code, self.response.text))

  def test_CorrectTemplate(self):
    '''Checks that the expected template is loaded'''
    for testCode in self.testCodes:
      self.assertIn(testCode, self.response.text, 'There appears to be a problem with the rendering of the template, TestCode: {0}. Template returns the following:\n{1}'.format(testCode, self.response.text))

  def test_ValidHtml(self):
    '''Checks HTML validity'''
    self.validate()
    responseData = json.loads(self.validationResponse.content)
    testPassed = True
    failureMessages = ''
    for message in responseData['messages']:
      if message['type'] == 'info':
        if 'subtype' in message.keys():
           if message['subtype'] == 'warning':
              testPassed = False
              failureMessages += self.constructFailureMessage(message)
      if message['type'] in ('error', 'non-document-error'):
        failureMessages += self.constructFailureMessage(message)
        testPassed = False
    self.assertTrue(testPassed, failureMessages)

class PostHttpTest(GetHttpTest):
  '''A test for post requests that do not use sessions'''

  _payload = {}
  '''The data to be POSTed to the sever'''

  def __init__(self, *args, **kwargs):
    super(PostHttptest, self).__init__(*args, **kwargs)
    self.payload = self._payload.copy()

  def setUp(self):
    self.response = requests.post(self.url, data=self.payload, params=self.params)

class GetHttpSessionTest(GetHttpTest):

  def __init__(self, *args, **kwargs):
    super(GetHttpSessionTest, self).__init__(*args, **kwargs)
    self.s = requests.Session()

  def setUp(self):
    self.response = self.s.get(self.url, params=self.params)

class PostHttpSessionTest(PostHttpTest):
  '''A test for post requests that use sessions (e.g. get decorated with logsInAs)'''

  def __init__(self, *args, **kwargs):
    super(PostHttpSessionTest, self).__init__(*args, **kwargs)
    self.s = requests.Session()

  def setUp(self):
    self.response = self.s.post(self.url, data=self._payload, params=self.params) 

class OneRedirectionMixin:
  '''A mixin for testing redirection pages.'''

  def test_redirect(self):
    '''Checks the response history for 302 redirects'''
    self.assertEqual(302, self.response.history[0].status_code)

def usesCsrf(c):
  '''A class decorator to indicate the test utilises csrf'''
  
  _oldSetup = c.setUp

  def setUp(self):
    getResponse = self.s.get(self.url)
    self.csrf = self.s.cookies.get_dict()['csrftoken']
    if hasattr(self, 'payload'):
      self.payload['csrfmiddlewaretoken'] = self.csrf #special case for post classes
    _oldSetup(self)

  c.setUp = setUp

  return c

def logsInAs(username, password, csrf=True):
  '''A class decorator that creates and logs in a user on setup, and deletes it on teardown. Should be applied BEFORE usesCsrf decorator'''

  def _logsInAs(c):
  
    c.loginUrl = c.baseUrl + reverse('login')
    _oldSetup = c.setUp
    _oldTearDown = c.tearDown
  
    def setUp(self):
      User.objects.create_user(username=username, password=password)
      if self.s is not None:
        self.s = requests.Session()
      getResponse = self.s.get(self.loginUrl)
      loginCsrf = self.s.cookies.get_dict()['csrftoken']
      loginResponse = self.s.post(self.loginUrl, data={'username':username, 'password':password, 'csrfmiddlewaretoken':loginCsrf})
      _oldSetUp(self)
  
    def tearDown(self):
      _oldTearDown(self)
      User.objects.filter(username=username).delete()

    c.setUp = setUp
    c.tearDown = tearDown
    return c
    
  return _logsInAs

def signsExampleLicense(username):
  '''A class decorator that creates a test license and makes the user specified by username sign it on setUp'''
  def _signsExampleLicense(c):

    _oldSetup = c.setUp
    _oldTearDown = c.tearDown
     
    license = License(text='This is an example license used in a test', effectiveDate=date.today() - timedelta(1))

    def setUp(self)
      user = User.objects.get(username=username)
      license.save()
      self.agreement = LicenseAgreement(user=User, text=license)
      _oldSetup(self)

    def tearDown(self):
      self.agreement.delete()
      license.delete()
      _oldTearDown(self)

    c.setUp = setUp
    c.tearDown = tearDown
    return c
  return _signsExampleLicense


def joinsLabGroup(username, labGroupTitle):
  '''A class decorator that creates a test lab group with labGroupTitle as it's title and assigns user identified by
  username to that lab group'''
  def _joinsLabGroup(c):
    _oldSetup = c.setUp
    _oldTearDown = c.tearDown

    labGroup = LabGroup(labGroupTitle, 'War drobe', 'Aslan@example.com', 'new_magic')

    def setUp(self):
      user = User.objects.get(username=username)
      user.labgroup_set.add(labGroup)
      labGroup.save()
      _oldSetup(self)

    def tearDown(self):
      _oldTearDown(self)
      labGroup.delete()

    c.setUp = setUp
    c.tearDown = tearDown
    return c
  return _joinsLabGroup

def createsChemicalClass(label, description):
  '''A class decorator that creates a test chemical class for the addition of compounds into the database'''

  def _createsChemicalClass(c):

    _oldSetup = c.setUp
    _oldTearDown = c.tearDown

    chemicalClass = ChemicalClass(label, description)

    def setUp(self):
      chemicalClass.save()
      _oldSetup(self)

    def tearDown(self):
      chemicalClass.delete()
      _oldTearDown(self)

    c.setUp = setUp
    c.tearDown = tearDown
    return c
  return _createsChemicalClass

def choosesLabGroup(username, labGroupTitle):
  '''A class decorator that sets up a user to be using a given labgroup for a session to view e.g. compound lists
  Necessarily assumes that the user has been logged in and adjoined to the lab group
  '''

  def _choosesLabGroup(username, labGroupTitle):

    _oldSetUp = c.setUp
    _oldTearDown = c.tearDown
    c.groupSelectUrl = c.baseUrl + reverse('selectGroup')

    def setUp(self):
      labGroup = LabGroup.objects.get(title=labGroupTitle)
      user = User.objects.get(username)
      getResponse = self.s.get(self.groupSelectUrl)
      selectCsrf = self.s.gookies.get_dict()['csrftoken']
      self.s.post(self.groupSelectUrl, data={'labGroup':labGroup.id, 'csrfmiddlewaretoken':selectCsrf})
      _oldSetUp(self)

    c.setUp = setUp
    return c
  return _choosesLabGroup
