#!/usr/bin/env python
"""This module contains tests for teh confirmation page.."""

from HttpTest import GetHttpTest, PostHttpTest, PostHttpSessionTest, usesCsrf
from django.contrib.auth.models import User
from DRP.models import ConfirmationCode
from uuid import uuid4
import requests
import unittest
from DRP.tests import runTests

loadTests = unittest.TestLoader().loadTestsFromTestCase


class ConfirmationPage(GetHttpTest):

    """Tests the simple GET request."""

    _params = {'code': uuid4()}
    url = GetHttpTest.baseUrl + '/confirm.html'
    testCodes = ['340545b4-9577-489c-b7ec-e664d8bbbe4c']


class ConfirmationPage2(GetHttpTest):

    """Tests the simple GET request without the required URL parameter (403)."""

    url = GetHttpTest.baseUrl + '/confirm.html'
    testCodes = ['b4da5e80-190b-4fe4-a97c-7f8bb9c213a5']

    def test_Status(self):
        """Ensure that response code correctly 403 when missing URL parameter."""
        self.assertEqual(403, self.response.status_code)


@usesCsrf
class PostConfirmationPage(PostHttpSessionTest):

    """Tests the POST request for a user with correct credentials in the correct state."""

    confirmationCode = uuid4()
    _params = {'code': confirmationCode}
    url = PostHttpTest.baseUrl + '/confirm.html'
    testCodes = ['9a301eb9-1f13-43ef-b3eb-3d66ebd4d566']

    def setUp(self):
        """Test the POST request for a user with correct credentials in the correct state."""
        self.user = User.objects.create_user(username='Aslan', password='banana', email='aslan@example.com')
        self.user.is_active = False
        self.user.save()
        self.code = ConfirmationCode(user=self.user, code=self.confirmationCode)
        self.code.save()
        self.response = self.s.post(self.url, data={'username': 'Aslan', 'password': 'banana', 'csrfmiddlewaretoken': self.csrf}, params={'code': self.code.code})

    def tearDown(self):
        """Delete objects created for this test."""
        self.user.delete()
        self.code.delete()


@usesCsrf
class PostConfirmationPage2(PostConfirmationPage):

    """Test for a 403 with non-matching confirmation code."""

    confirmationCode = uuid4()
    testCodes = ['b4da5e80-190b-4fe4-a97c-7f8bb9c213a5']

    def setUp(self):
        """Set up active use to look for POST confirmation response."""
        self.user = User.objects.create_user(username='Aslan', password='banana', email='aslan@example.com')
        self.user.is_active = False
        self.user.save()
        self.code = ConfirmationCode(user=self.user, code=self.confirmationCode)
        self.code.save()
        self.response = self.s.post(self.url, data={'username': self.user.username, 'password': 'banana', 'csrfmiddlewaretoken': self.csrf}, params={'code': uuid4()})

    def test_Status(self):
        """Look for POST confirmation response."""
        self.assertEqual(403, self.response.status_code)

suite = unittest.TestSuite([
    loadTests(ConfirmationPage),
    loadTests(ConfirmationPage2),
    loadTests(PostConfirmationPage),
    loadTests(PostConfirmationPage2)
])

if __name__ == '__main__':
    runTests(suite)
