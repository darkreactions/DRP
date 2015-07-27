#!/usr/bin/env python
'''This file contains tests which validate the html output and will eventually
also contain unit tests for javascript components.

These tests should NOT be utilised for testing django functionality of webpages, which
should be managed by other unit tests on views (worst case) or their components
(preferred).
'''

import unittest
import TestConfig
import requests
import json
from DRP.settings import TESTING_SERVER_NAME, EXTERNAL_HTML_VALIDATOR
loadTests = unittest.TestLoader.loadTestsFromTestCase

class HomePage(unittest.TestCase):
  #This class exemplifies the standard structure of a test. Check the documentation for 'rolling your own'

  baseUrl = 'http://' + TESTING_SERVER_NAME

  def setUp(self):
    self.response = requests.get(self.baseUrl)

  def validate(self):
    self.validationResponse = requests.post(EXTERNAL_HTML_VALIDATOR, params={'out':'json'}, data=self.response.content, headers={'content-type':self.response.headers.get('content-type')})

  @staticmethod
  def constructFailureMessage(message):
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
    self.assertEqual(self.response.status_code, 200)

  def test_ValidHtml(self):
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

class AboutPage(HomePage):

  def setUp(self):
    self.response = requests.get(self.baseUrl + '/about.html')
  
def suite():
  #This function should be adjusted to contain the loadTests() function enacted on each test case.
  return unittest.TestSuite([
          loadTests(HomePage),
          loadTests(AboutPage)
          ])

if __name__ == '__main__':
  #Runs the test- a good way to check that this particular test set works without having to run all the tests.
  unittest.main()
