#!/usr/bin/env python
'''This contains tests for the registration page'''

import time
import unittest
import requests
from django.contrib.auth.models import User
from DRP.models import ConfirmationCode
from ContactPage import ContactPage, ContactPage_POST, usesCsrf
from django.conf import settings
from django.core.urlresolvers import reverse
import imaplib
import email
import re
import string
import random

loadTests = unittest.TestLoader().loadTestsFromTestCase
#The registration page tests assume that the django provided form will behave correctly.

class RegisterPage(ContactPage):
  '''Checks the register page html validity'''

  url=ContactPage.baseUrl + '/register.html'
  testCodes = ['4cf1abe0-9118-471c-ac4e-34e863e87402','be088572-3adc-4757-8059-d16db2ea77a6'] 

class RegisterPage_POST(ContactPage_POST):
  '''Checks the register page email response'''
    
  url = ContactPage_POST.baseUrl + '/register.html'
  templateId = 'd776703c-bf1c-4a0a-89d1-1fcd83093967'
  emailCode = 'ff5327c2-9e62-4f4b-b64c-67bd9e9935c2'
  emailHeader = 'Dark Reactions Project Registration'
  testCodes = ['8f7aa6e8-2be5-4630-a205-5fceeea81fe4']

  @usesCsrf
  def setUp(self):
    self.username = ''.join(random.choice(string.ascii_uppercase + string.ascii_lowercase) for i in range(6))
    self.response = self.s.post(self.url, data={'username':'testing', 'first_name':self.username, 'last_name':'user', 'password1':'testpass', 'password2':'testpass', 'email':settings.EMAIL_HOST_USER, 'csrfmiddlewaretoken':self.csrf})

  def tearDown(self):
    users = User.objects.filter(username='testing')
    for user in users:
      user.delete()

  def emailCheck(self, emailText):
    '''Worker method for test_email, actually does the checks, but allows test_email to do the heavy lifting of actually fetching the mail'''
    linkRe = re.compile('http://' + settings.SERVER_NAME + reverse('confirm') + '\?code=(?P<code>\w{8}-\w{4}-\w{4}-\w{4}-\w{12})')
    match = linkRe.search(emailText) 
    if self.emailCode in emailText and self.username in emailText:
      if match: 
        try:
          self.assertIsNotNone(match, 'The link could not be found in the email')
          self.code = ConfirmationCode.objects.get(code=match.group('code'))
          self.assertEqual(settings.EMAIL_HOST_USER, self.code.user.email, 'the email address is somehow wrong...')
          user = self.code.user
          r = requests.get('http://' + settings.SERVER_NAME + reverse('confirm') + '?code={0}'.format(match.group('code')))
          self.assertEqual(r.status_code, 200, 'The link from the email returned a wrong status code')
          self.assertIsNotNone(user)
          self.assertFalse(user.is_active)
        except ConfirmationCode.DoesNotExist as e:
          codes = ', '.join([c.code for c in ConfirmationCode.objects.all()])
          self.assertTrue(False, 'Could not find code. Codes available: {0}, This code: {1}'.format(codes, match.group('code')))
      else:
        self.assertTrue(False, '{0} did not match {1}'.format(linkRe.pattern, emailText))
      return True
    else:
      return False
  
  def test_email(self):
      '''Checks that the registration email has been recieved and is in good shape.''' 
      messages = ''
      testPass = False
      time.sleep(10)
      m = imaplib.IMAP4_SSL(settings.EMAIL_IMAP_HOST)
      m.login(settings.EMAIL_HOST_USER, settings.EMAIL_HOST_PASSWORD)
      sel_rv, data = m.select(settings.EMAIL_IMAP_INBOX)
      try:
        if sel_rv=='OK':
          search_rv, data = m.search(None, 'FROM', settings.DEFAULT_FROM_EMAIL, 'SUBJECT', '{0}'.format(self.emailHeader))
          if search_rv=='OK':
            if len(data[0].split()) > 0:
              for num in data[0].split():
                fetch_rv, msgData = m.fetch(num, '(RFC822)')
                if fetch_rv=='OK':
                    if self.emailCheck(str(email.message_from_string(msgData[0][1]))): #This check should do something to make sure that the email is unique
                      m.store(num, '+FLAGS','\\DELETED')
                      m.expunge()
                else:
                  messages += str(email.message_from_string(msgData[0][1]))
            else:
              raise RuntimeError('Message not found in inbox')
          else:
            raise RuntimeError('Message Searching failed')
        else:
          raise RuntimeError('Inbox selection failed. Perhaps a different inbox is needed for settings.EMAIL_IMAP_INBOX in settings.py.')
      finally:
        m.close()
        m.logout()

tests = [loadTests(RegisterPage)]

if settings.TESTING:
  tests.append(loadTests(RegisterPage_POST))

suite = unittest.TestSuite(tests)

if __name__ == '__main__':
  #Runs the test- a good way to check that this particular test set works without having to run all the tests.
  unittest.main()
