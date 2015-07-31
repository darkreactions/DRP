#!/usr/bin/env python
'''This contains tests for the registration page'''

import unittest
import requests
from django.contrib.auth.models import User
from LoginPage import LoginPage
from ContactPage import usesCsrf

#The registration page tests assume that the django provided form will behave correctly.

class RegisterPage(LoginPage):
  '''Checks the register page html validity'''

  url= LoginPage.baseUrl + '/register.html'

class RegisterPage_POST(LoginPage_POST):
  '''Checks the register page POST redirect for redirecting to the login page'''
    
  url = LoginPage_POST.baseUrl + '/register.html'
  templateId = 'd776703c-bf1c-4a0a-89d1-1fcd83093967'

  @usesCsrf
  def setUp(self):
    self.email = 'aslan@example.com'
    self.response = self.s.post(self.url, data={'username':"testUser", 'password1':'testpass', 'password2':'testpass', 'email':self.email, 'csrfmiddlewaretoken':self.csrf})

  def tearDown(self):
    users = User.objects.filter(email=self.email)
    for user in users:
      user.delete() 
    
def Suite():
  return unittest.TestSuite([
   loadtests(AboutPage),
   loadtests(ContactPage),
   loadtests(ContactPage_POST),
   loadtests(RegisterPage),
   loadtests(RegisterPage_POST),
  ])

if __name__ == '__main__':
  #Runs the test- a good way to check that this particular test set works without having to run all the tests.
  unittest.main()
