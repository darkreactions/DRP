#!/usr/bin/env python
'''Tests for the Email classes contained in DRP.Email'''

import TestConfig
import unittest
from DRP.Email import Email
from DRP.settings import EMAIL_HOST_USER, EMAIL_HOST_PASSWORD, EMAIL_IMAP_HOST, DEFAULT_FROM_EMAIL, SKIP_EMAIL_TESTS, EMAIL_IMAP_INBOX 
from uuid import uuid4
import imaplib
import time
import email

loadTests = unittest.TestLoader.loadTestsFromTestCase

class EmailSendsAndRecieves(unittest.TestCase):
  '''Sends and recieves a test email, and checks that the contents are correct'''

  def setUp(self):
      '''Create an email with a unique header'''
      self.headerId = uuid4()
      self.email = Email('Test Subject Header: {0}'.format(self.headerId), 'This message is a test. Please disregard but do not delete this email.', [EMAIL_HOST_USER])

  def runTest(self):
      '''Sends an email using SMTP and fetches it via IMAP'''
      errMessage = 'False is not True'
      messages = ''
      testPass = False
      self.email.send()
      time.sleep(10)
      m = imaplib.IMAP4_SSL(EMAIL_IMAP_HOST)
      m.login(EMAIL_HOST_USER, EMAIL_HOST_PASSWORD)
      sel_rv, data = m.select(EMAIL_IMAP_INBOX)
      try:
        if sel_rv=='OK':
          search_rv, data = m.search(None, 'FROM', DEFAULT_FROM_EMAIL, 'SUBJECT', 'Test Subject Header: {0}'.format(self.headerId))
          if search_rv=='OK':
            if len(data[0].split()) > 1:
              raise RuntimeError('More than one message with unique ID. Test aborted.')
            else:
              for num in data[0].split():
                fetch_rv, msgData = m.fetch(num, '(RFC822)')
                if fetch_rv=='OK':
                    testPass = True
                    m.store(num, '+FLAGS','\\DELETED')
                    m.expunge()
                else:
                  messages += str(email.message_from_string(msgData[0][1]))
              else:
                raise RuntimeError('Message Fetching failed')
          else:
            raise RuntimeError('Message Searching failed')
        else:
          raise RuntimeError('Inbox selection failed. Perhaps a different inbox is needed for EMAIL_IMAP_INBOX in settings.py.')
      except RuntimeError as e:
        errMessage = repr(e)
      finally:
        m.close()
        m.logout()
      self.assertTrue(testPass, errMessage + messages)


def suite():
  if not SKIP_EMAIL_TESTS:
    return unittest.TestSuite([
          loadTests(EmailSendsAndRecieves)
          ])
  else:
    return unittest.TestSuit([])

if __name__ == '__main__':
  #Runs the test- a good way to check that this particular test set works without having to run all the tests.
  unittest.main()
