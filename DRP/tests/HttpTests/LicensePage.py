#!/usr/bin/env python
'''This module provides tests for the License agreement page'''

import unittest
from django.contrib.auth.models import User
from DRP.models import License, LicenseAgreement
import requests
import datetime
from HttpTest import GetHttpTest, PostHttpTest, GetHttpSessionTest, PostHttpSessionTest, usesCsrf, logsInAs
loadTests = unittest.TestLoader().loadTestsFromTestCase


@logsInAs('Aslan', 'banana')
class LicenseAgreementPage(GetHttpSessionTest):

    url = GetHttpTest.baseUrl + '/license.html'
    testCodes = ['c9e46ba1-cd2a-4080-88b5-97415fa7c484']

    def setUp(self):
        '''Sets up the test by requesting the page uri'''
        self.license = License(text='some test', effectiveDate=datetime.date.today() - datetime.timedelta(1))
        self.license.save()
        self.response = self.s.get(self.url, params=self.params)

    def tearDown(self):
        self.license.delete()


@logsInAs('Aslan', 'banana')
@usesCsrf
class PostLicenseAgreementPage(PostHttpSessionTest):
    '''defines a test case with good credentials leading to a redirect'''

    url = PostHttpTest.baseUrl + '/license.html'
    testCodes = ['3a9f74ee-5c78-4ec0-8893-ce0476808131']

    def setUp(self):
        self.license = License(text='some test', effectiveDate=datetime.date.today() - datetime.timedelta(1))
        self.license.save()
        self.response = self.s.post(self.url, data={'username': 'Aslan', 'password': 'banana', 'licenseId': self.license.id, 'csrfmiddlewaretoken': self.csrf}, params={'next': '/contact.html'})

    def test_Redirect(self):
        self.assertEqual(302, self.response.history[0].status_code)

    def tearDown(self):
        self.license.delete()


@logsInAs('Aslan', 'banana')
@usesCsrf
class PostLicenseAgreementPage2(PostHttpSessionTest):
    '''Defines a test case for good credentials with no redirect'''

    url = PostHttpTest.baseUrl + '/license.html'
    testCodes = ['9d6147f1-e321-4aff-8d04-966ca24a2ab0']

    def setUp(self):
        self.license = License(text='some test', effectiveDate=datetime.date.today() - datetime.timedelta(1))
        self.license.save()
        self.response = self.s.post(self.url, data={'username': 'Aslan', 'password': 'banana', 'licenseId': self.license.id, 'csrfmiddlewaretoken': self.csrf})

    def tearDown(self):
        self.license.delete()


@logsInAs('Aslan', 'banana')
@usesCsrf
class PostLicenseAgreementPage3(PostLicenseAgreementPage):
    '''Defines a testcase where the licenseagreement is already signed, with a redirect'''

    def setUp(self):
        self.user = User.objects.get(username='Aslan')
        self.license = License(text='some test', effectiveDate=datetime.date.today() - datetime.timedelta(1))
        self.license.save()
        self.agreement = LicenseAgreement(user=self.user, text=self.license)
        self.agreement.save()
        self.response = self.s.post(self.url, data={'username': 'Aslan', 'password': 'banana', 'csrfmiddlewaretoken': self.csrf}, params={'next': '/contact.html'})

    def tearDown(self):
        self.agreement.delete()
        self.license.delete()


@logsInAs('Aslan', 'banana')
@usesCsrf
class PostLicenseAgreementPage4(PostLicenseAgreementPage2):
    '''Defines a testcase where the licenseagreement is already signed, without a redirect'''

    testCodes = ['c87f2095-c8b2-4798-89bd-2b93ee33d338']

    def setUp(self):
        self.user = User.objects.get(username='Aslan')
        self.license = License(text='some test', effectiveDate=datetime.date.today() - datetime.timedelta(1))
        self.license.save()
        self.agreement = LicenseAgreement(user=self.user, text=self.license)
        self.agreement.save()
        self.response = self.s.post(self.url, data={'username': 'Aslan', 'password': 'banana', 'csrfmiddlewaretoken': self.csrf})

    def tearDown(self):
        self.agreement.delete()
        self.license.delete()

suite = unittest.TestSuite([
    loadTests(LicenseAgreementPage),
    loadTests(PostLicenseAgreementPage),
    loadTests(PostLicenseAgreementPage2),
    loadTests(PostLicenseAgreementPage3),
    loadTests(PostLicenseAgreementPage4),
])

if __name__ == '__main__':
    unittest.TextTestRunner(verbosity=2).run(suite)
