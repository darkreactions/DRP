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
import html5lib
from DRP.settings import TESTING_SERVER_NAME, EXTERNAL_HTML_VALIDATOR, EMAIL_HOST_USER
loadTests = unittest.TestLoader.loadTestsFromTestCase

class HomePage(unittest.TestCase):
  '''Tests the home page for HTML validity'''

  baseUrl = 'http://' + TESTING_SERVER_NAME
  url = baseUrl

  def setUp(self):
    '''Sets up the test by requesting the home page uri'''
    self.response = requests.get(self.baseUrl)

  def validate(self):
    '''Sends the output of the requested page to the w3c html validator to validate'''
    self.validationResponse = requests.post(EXTERNAL_HTML_VALIDATOR, params={'out':'json'}, data=self.response.content, headers={'content-type':self.response.headers.get('content-type')})

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
    self.assertEqual(self.response.status_code, 200, 'Url {0} returns code {1}'.format(self.url, self.response.status_code))

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

class AboutPage(HomePage):
  '''Performs a get request on the AboutPage to check for Html Validity'''

  url = HomePage.baseUrl + '/about.html'

  def setUp(self):
    self.response = requests.get(self.url)

class ContactPage(AboutPage):
  '''Performs a simple get request on the contact page to test html validity'''

  url = AboutPage.baseUrl + '/contact.html'

class ContactPage_POST(ContactPage):
  '''Performs a correctly formed POST request to the contact page to test html validity'''
  
  def setUpCsrf(self):
    '''Obtains the relevant CSRF token from the page'''
    self.s = requests.Session()
    getResponse = self.s.get(self.url)
    self.csrf = self.s.cookies.get_dict()['csrftoken']

  def setUp(self):
    self.setUpCsrf()
    self.response = self.s.post(self.url, data={'email':EMAIL_HOST_USER, 'content':'This is a test message.', 'csrfmiddlewaretoken':self.csrf})

class ContactPage_POST_bad(ContactPage_POST):
  '''Perfoems a badly formed POST request to the contact page to test html validity'''

  def setUp(self):
    self.setUpCsrf()
    self.response = self.s.post(self.url, data={'email':EMAIL_HOST_USER, 'content':'', 'csrfmiddlewaretoken':self.csrf})
  
class ContactPage_POST_bad2(ContactPage_POST):
  '''Performs a badly formed POST request missing a field entirely to test html validity'''

  def setUp(self):
    self.setUpCsrf()
    self.response = self.s.post(self.url, data={'content':'This is a test.', 'csrfmiddlewaretoken':self.csrf})

class ContactPage_POST_bad3(ContactPage):
  '''Confirms that CSRF has been correctly implemented by forgetting the csrf token.'''

  def setUp(self):
    self.response = requests.post(self.url, data={'content':'This is a test.'})

  def test_Status(self):
    self.assertEqual(403, self.response.status_code)

def suite():
  #This function should be adjusted to contain the loadTests() function enacted on each test case.
  return unittest.TestSuite([
          loadTests(HomePage),
          loadTests(AboutPage),
          loadTests(ContactPage),
          loadTests(ContactPage_POST),
          loadTests(ContactPage_POST_bad),
          loadTests(ContactPage_POST_bad2),
          loadTests(ContactPage_POST_bad3)
          ])

if __name__ == '__main__':
  #Runs the test- a good way to check that this particular test set works without having to run all the tests.
  unittest.main()
