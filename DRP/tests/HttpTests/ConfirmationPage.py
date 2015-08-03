#!/usr/bin/env python
'''This module contains tests for teh confirmation page'''

from ContactPage import ContactPage, ContactPage_POST, usesCsrf
from django.contrib.auth.models import User
from DRP.models import ConfirmationCode
from uuid import uuid4
import requests
import unittest

loadTests = unittest.TestLoader().loadTestsFromTestCase

class ConfirmationPage(ContactPage):
  '''Tests the simple GET request'''

  url = ContactPage.baseUrl + '/confirm.html?code={0}'.format(uuid4())
  testCodes = ['340545b4-9577-489c-b7ec-e664d8bbbe4c']

class ConfirmationPage2(ContactPage):
  '''Tests the simple GET request without the required URL parameter (403)'''

  url = ContactPage.baseUrl + '/confirm.html?'
  testCodes = ['b4da5e80-190b-4fe4-a97c-7f8bb9c213a5']
  
  def test_Status(self):
    self.assertEqual(403, self.response.status_code)

class ConfirmationPage_POST(ContactPage_POST):
  '''Tests the POST request for a user with correct credentials in the correct state'''
 
  confirmationCode =uuid4() 
  url = ContactPage.baseUrl + '/confirm.html?code={0}'.format(confirmationCode)
  testCodes = ['9a301eb9-1f13-43ef-b3eb-3d66ebd4d566']

  @usesCsrf
  def setUp(self):
    self.user = User.objects.create_user(username='Aslan', password='banana', email='aslan@example.com')
    self.user.is_active = False
    self.user.save()
    self.code = ConfirmationCode(user=self.user, code=self.confirmationCode)
    self.code.save()
    self.response = self.s.post(self.url, data={'username':'Aslan', 'password':'banana', 'csrfmiddlewaretoken':self.csrf}, params={'code':self.code.code})
  
  def tearDown(self):
    self.user.delete()
    self.code.delete()


class ConfirmationPage_POST2(ConfirmationPage_POST):
  '''Tests for a 403 with non-matching confirmation code'''

  confirmationCode = uuid4()
  testCodes = ['b4da5e80-190b-4fe4-a97c-7f8bb9c213a5']

  @usesCsrf
  def setUp(self):
    self.user = User.objects.create_user(username='Aslan', password='banana', email='aslan@example.com')
    self.user.is_active = False
    self.user.save()
    self.code = ConfirmationCode(user=self.user, code=self.confirmationCode)
    self.code.save()
    self.response = self.s.post(self.url + '?code={0}'.format(uuid4()), data={'username':self.user.username, 'password':'banana', 'csrfmiddlewaretoken':self.csrf})

  def test_Status(self):
    self.assertEqual(403, self.response.status_code)

suite = unittest.TestSuite([
  loadTests(ConfirmationPage),
  loadTests(ConfirmationPage2),
  loadTests(ConfirmationPage_POST),
  loadTests(ConfirmationPage_POST2)
])

if __name__ == '__main__':
  unittest.main()
