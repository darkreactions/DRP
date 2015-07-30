#!/usr/bin/env python
'''This package provides tests for the Contact page'''
#TODO: test for email sending

from AboutPage import AboutPage
import requests
import unittest
loadTests = unittest.TestLoader().loadTestsFromTestCase

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
    self.response = self.s.post(self.url, data={'email':'aslan@example.com', 'content':'This is a test message.', 'csrfmiddlewaretoken':self.csrf})

class ContactPage_POST_bad(ContactPage_POST):
  '''Perfoems a badly formed POST request to the contact page to test html validity'''

  def setUp(self):
    self.setUpCsrf()
    self.response = self.s.post(self.url, data={'email':'aslan@example.com', 'content':'', 'csrfmiddlewaretoken':self.csrf})
  
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

suite = unittest.TestSuite([
  loadTests(ContactPage),
  loadTests(ContactPage_POST_bad),
  loadTests(ContactPage_POST_bad2),
  loadTests(ContactPage_POST_bad3)
])

if __name__ == '__main__':
  #Runs the test- a good way to check that this particular test set works without having to run all the tests.
  unittest.main()
