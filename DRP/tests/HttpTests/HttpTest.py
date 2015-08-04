'''A module containing base classes and decorators for HttpTests'''

import requests
from django.conf import settings
from django.core.urlresolvers import reverse
from django.contrib.auth.models import User
from abc import ABCMeta
from DRP.tests import DRPTestCase
import json

class GetHttpTest(DRPTestCase):

  baseUrl = 'http://' + settings.SERVER_NAME
  url = baseUrl
  testCodes = []
  params = {}
  '''any GET params to be added to the reuqest.'''

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
    self.assertEqual(self.response.status_code, 200, 'Url {0} returns code {1}. Page content follows:\n\n{2}'.format(self.url, self.response.status_code, self.response.text))

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

  payload = {}
  '''The data to be POSTed to the sever'''

  def setUp(self):
    self.response = requests.post(self.url, data=self.payload, params=self.params)

def usesCsrf(c):
  '''A class decorator to indicate the test utilises csrf'''
  class CsrfWrapped(c):

    s=None

    def setUp(self):
      if not self.s:
        self.s = requests.Session()
      getResponse = self.s.get(self.url)
      self.csrf = self.s.cookies.get_dict()['csrftoken']
      super(CsrfWrapped, self).setUp()

  return CsrfWrapped

def logsInAs(username, password, csrf=True):
  '''A class decorator that creates and logs in a user on setup, and deletes it on teardown. Should be applied BEFORE usesCsrf decorator'''

  def _logsInAs(c):
    class LogsInWrapped(c):
  
      loginUrl = c.baseUrl + reverse('login')
      s = None 
  
      def setUp(self):
        User.objects.create_user(username=username, password=password)
        if not self.s:
          self.s = requests.Session()
        getResponse = self.s.get(self.loginUrl)
        loginCsrf = self.s.cookies.get_dict()['csrftoken']
        loginResponse = self.s.post(self.loginUrl, data={'username':username, 'password':password, 'csrfmiddlewaretoken':loginCsrf})
        super(LogsInWrapped, self).setUp()
  
      def tearDown(self):
        super(LogsInWrapped, self).tearDown()
        User.objects.filter(username=username).delete()
    
    return LogsInWrapped
  return _logsInAs
