#/usr/bin/env python
'''Provides HTML tests for the login'''
import unittest
import requests
from ContactPage import ContactPage_POST, usesCsrf
from AboutPage import AboutPage

loadTests = unittest.TestLoader().loadTestsFromTestCase

#NOTE the correct behaviour of the default login view by django response is taken for granted.

class LoginPage(AboutPage):
  '''Confirms that GETing the login page results in valid html'''

  url = AboutPage.baseUrl + '/login.html'

  def setUp(self):
    self.response = requests.get(self.url)

class LoginPage_POST(ContactPage_POST):
  '''confirms that posting valid data to the login page results in a redirect''' 

  url = ContactPage_POST.baseUrl + '/login.html'
  testCodes = [u'3a9f74ee-5c78-4ec0-8893-ce0476808131', '1f47e7ab-1900-4683-ba1c-63330ec2f71a']

  @usesCsrf
  def setUp(self):
    self.tmpUser = User.objects.create_user(username="testUser", email=settings.EMAIL_HOST_USER, password="testpass")
    self.tmpUser.save()
    self.response = self.s.post(self.url, data={'username':"testUser", 'password':"testpass", 'csrfmiddlewaretoken':self.csrf}, params={'next':'/contact.html'})

  def test_Status(self):
    self.assertEqual(200, self.response.status_code)
    self.assertEqual(302, self.response.history[0].status_code)

  def tearDown(self):
    self.tmpUser.delete()

suite = unittest.TestSuite([
  loadtests(LoginPage_POST),
  loadtests(LoginPage),
])

if __name__ == '__main__':
  #Runs the test- a good way to check that this particular test set works without having to run all the tests.
  unittest.main()
