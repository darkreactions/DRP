#!/usr/bin/env python
"""Provides HTML tests for the login."""

import unittest
import requests
from HttpTest import GetHttpTest, PostHttpTest, usesCsrf, PostHttpSessionTest
from django.contrib.auth.models import User
from django.conf import settings
from DRP.tests import runTests

loadTests = unittest.TestLoader().loadTestsFromTestCase

# NOTE the correct behaviour of the default login view by django response
# is taken for granted.


class LoginPage(GetHttpTest):

    """Confirms that GETing the login page results in valid html."""

    url = GetHttpTest.baseUrl + '/login.html'


@usesCsrf
class PostLoginPage(PostHttpSessionTest):

    """Confirms that posting valid data to the login page results in a redirect."""

    url = PostHttpTest.baseUrl + '/login.html'
    testCodes = [u'3a9f74ee-5c78-4ec0-8893-ce0476808131',
                 '1f47e7ab-1900-4683-ba1c-63330ec2f71a']

    def setUp(self):
        """Make POST request of valid data to be logged in."""
        self.tmpUser = User.objects.create_user(
            username="testUser", email=settings.EMAIL_HOST_USER, password="testpass")
        self.tmpUser.save()
        self.response = self.s.post(self.url, data={
                                    'username': "testUser", 'password': "testpass", 'csrfmiddlewaretoken': self.csrf}, params={'next': '/contact.html'})

    def test_Status(self):
        """Ensure that the test redirects to the correct response."""
        self.assertEqual(200, self.response.status_code)
        self.assertTrue(len(self.response.history) > 0, self.response.content)
        self.assertEqual(302, self.response.history[
                         0].status_code, self.response.content)

    def tearDown(self):
        """Delete the temporary user for this test."""
        self.tmpUser.delete()

suite = unittest.TestSuite([
    loadTests(LoginPage),
    loadTests(PostLoginPage),
])

if __name__ == '__main__':
    # Runs the test- a good way to check that this particular test set works
    # without having to run all the tests.
    runTests(suite)
