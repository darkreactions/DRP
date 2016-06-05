#!/usr/bin/env python
'''This module provides tests for the Lab group agreement page'''

import unittest
from django.contrib.auth.models import User
from DRP.models import LabGroup
import requests
from django.core.urlresolvers import reverse
from HttpTest import GetHttpTest, PostHttpTest, GetHttpSessionTest, PostHttpSessionTest, usesCsrf, logsInAs
from DRP.tests import runTests
from DRP.tests.decorators import signsExampleLicense
loadTests = unittest.TestLoader().loadTestsFromTestCase


@logsInAs('Aslan', 'banana')
@signsExampleLicense('Aslan')
class LabGroupPage(GetHttpSessionTest):
    '''Gets the page when we have lab groups'''

    url = GetHttpTest.baseUrl + reverse('joinGroup')
    testCodes = ['660d4702-1621-406a-a78e-fe34c8a5a721']

    def setUp(self):
        '''Sets up the test by requesting the page uri'''
        self.labGroup = LabGroup.objects.makeLabGroup('test', 'narnia', 'Aslan@example.com', 'old_magic')
        self.labGroup.save()
        self.response = self.s.get(self.url, params=self.params)

    def tearDown(self):
        self.labGroup.delete()


@logsInAs('Aslan', 'banana')
@signsExampleLicense('Aslan')
class LabGroupPage404(GetHttpSessionTest):
    '''Gets the page when we have no lab groups'''

    url = GetHttpTest.baseUrl + reverse('joinGroup')
    testCodes = ['660d4702-1621-406a-a78e-fe34c8a5a721']

    def setUp(self):
        '''Sets up the test by requesting the page uri'''
        self.response = self.s.get(self.url, params=self.params)

    def test_Status(self):
        self.assertEqual(self.response.status_code, 404, 'Url {0} returns code {1}. Page content follows:\n\n{2}'.format(self.url, self.response.status_code, self.response.text))


@logsInAs('Aslan', 'banana')
@signsExampleLicense('Aslan')
@usesCsrf
class PostLabGroupPage(PostHttpSessionTest):
    '''Posts to the page and checks that we joined the group'''

    url = PostHttpTest.baseUrl + reverse('joinGroup')
    testCodes = ['00c9fe70-a51f-4c9a-9d99-b88292ece120']

    def setUp(self):
        self.labGroup = LabGroup.objects.makeLabGroup('test', 'narnia', 'Aslan@example.com', 'old_magic')
        self.labGroup.save()
        self.response = self.s.post(self.url, data={'labGroup': self.labGroup.id, 'accessCode': 'old_magic', 'csrfmiddlewaretoken': self.csrf})

    def tearDown(self):
        self.labGroup.delete()

suite = unittest.TestSuite([
    loadTests(LabGroupPage),
    loadTests(PostLabGroupPage)
])

if __name__ == '__main__':
    runTests(suite)
