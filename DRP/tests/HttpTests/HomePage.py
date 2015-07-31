#!/usr/bin/env python
'''This package provides tests for the home page'''
from DRP.tests import DRPTestCase
from django.conf import settings
import unittest
import requests
import json

loadTests = unittest.TestLoader().loadTestsFromTestCase

class HomePage(DRPTestCase):
  '''Tests the home page for HTML validityi, and that the correct template is rendered.'''

  baseUrl = 'http://' + settings.SERVER_NAME
  url = baseUrl
  testCodes = ['2240b1ff-895c-458f-adf2-d04d85a164d1', 'd68a82db-bd18-4a9f-a1a2-03b3bb259595']
  '''A list of unique IDs. A check is performed to see that these are in the returned HTML'''

  def setUp(self):
    '''Sets up the test by requesting the home page uri'''
    self.response = requests.get(self.baseUrl)

  def validate(self):
    '''Sends the output of the requested page to the w3c html validator to validate'''
    self.validationResponse = requests.post(settings.EXTERNAL_HTML_VALIDATOR, params={'out':'json'}, data=self.response.content, headers={'content-type':self.response.headers.get('content-type')})

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

  def test_Status(self):
    '''Checks that the http response code is the expected value'''
    self.assertEqual(self.response.status_code, 200, 'Url {0} returns code {1}. Page content follows:\n\n{2}'.format(self.url, self.response.status_code, self.response.text))

  def test_CorrectTemplate(self):
    '''Checks that the expected template is loaded'''
    for testCode in self.testCodes:
      self.assertIn(testCode, self.response.text, 'There appears to be a problem with the rendering of the template, TestCode: {0}'.format(testCode))

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

suite = unittest.TestSuite([
  loadTests(HomePage)
])

if __name__ == '__main__':
  #Runs the test- a good way to check that this particular test set works without having to run all the tests.
  unittest.main()
