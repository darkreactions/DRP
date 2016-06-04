#!/usr/bin/env python
'''Tests for the Email classes contained in DRP.Email'''

from django.conf import settings
import unittest
from uuid import uuid4
import imaplib
import time
import email
import DRP.Email
from DRPTestCase import DRPTestCase, runTests

loadTests = unittest.TestLoader().loadTestsFromTestCase


class EmailSendsAndRecieves(DRPTestCase):
    '''Sends and recieves a test email, and checks that the contents are correct'''

    def setUp(self):
        '''Create an email with a unique header'''
        self.headerId = uuid4()
        self.email = DRP.Email.Email('Test Subject Header: {0}'.format(self.headerId), 'This message is a test. Please disregard but do not delete this email.', [settings.EMAIL_HOST_USER])

    def runTest(self):
        '''Sends an email using SMTP and fetches it via IMAP'''
        errMessage = 'False is not True'
        messages = ''
        testPass = False
        self.email.send()
        time.sleep(10)
        m = imaplib.IMAP4_SSL(settings.EMAIL_IMAP_HOST)
        m.login(settings.EMAIL_HOST_USER, settings.EMAIL_HOST_PASSWORD)
        sel_rv, data = m.select(settings.EMAIL_IMAP_INBOX)
        try:
            if sel_rv == 'OK':
                search_rv, data = m.search(None, 'FROM', settings.DEFAULT_FROM_EMAIL, 'SUBJECT', 'Test Subject Header: {0}'.format(self.headerId))
                if search_rv == 'OK':
                    for num in data[0].split():
                        fetch_rv, msgData = m.fetch(num, '(RFC822)')
                        if fetch_rv == 'OK':
                            testPass = True
                            m.store(num, '+FLAGS', '\\DELETED')
                            m.expunge()
                        else:
                            messages += str(email.message_from_string(msgData[0][1]))
                else:
                    raise RuntimeError('Message Searching failed')
            else:
                raise RuntimeError('Inbox selection failed. Perhaps a different inbox is needed for settings.EMAIL_IMAP_INBOX in settings.py.')
        except RuntimeError as e:
            errMessage = repr(e)
        finally:
            m.close()
            m.logout()
        self.assertTrue(testPass, errMessage + messages)

if not settings.SKIP_EMAIL_TESTS:
    suite = unittest.TestSuite([
        loadTests(EmailSendsAndRecieves)
    ])
else:
    suite = unittest.TestSuite([])


if __name__ == '__main__':
    # Runs the test- a good way to check that this particular test set works without having to run all the tests.
    runTests(suite)
