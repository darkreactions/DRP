#!/usr/bin/env python
"""This module provides tests for the License agreement page."""

import unittest
from django.contrib.auth.models import User
from DRP.models import LabGroup
from django.core.urlresolvers import reverse
import requests
from .httpTest import GetHttpTest, GetHttpSessionTest, PostHttpTest, logsInAs
loadTests = unittest.TestLoader().loadTestsFromTestCase


@logsInAs('Aslan', 'banana')
class AccountPageNoGroups(GetHttpSessionTest):
    """To test that the account page displays correctly with no group memberships."""

    url = GetHttpTest.baseUrl + reverse('account')
    testCodes = ['308fd1a4-d2da-4d4a-9a2b-58577e348050']

    def setUp(self):
        """Ensure that no group memberships visible."""
        self.response = self.s.get(self.url, params=self.params)


@logsInAs('Aslan', 'banana')
class AccountPageGroups(AccountPageNoGroups):
    """To test that the account page displays correctly when the user is a member of at least one group."""

    testCodes = ['48295bf1-5be1-4f94-aab6-1e5b7e97681b']

    def setUp(self):
        """Add user to lab group and then ensure that membership displayed properly."""
        self.labGroup = LabGroup.objects.makeLabGroup(
            'test', 'War Drobe', 'aslan@example.com', 'ancient_magic')
        self.labGroup.save()
        self.labGroup.users.add(User.objects.get(username='Aslan'))
        self.labGroup.save()
        self.response = self.s.get(self.url, params=self.params)

    def tearDown(self):
        """Remove created test lab group."""
        self.labGroup.delete()

suite = unittest.TestSuite([
    loadTests(AccountPageNoGroups),
    loadTests(AccountPageGroups)
])

if __name__ == '__main__':
    unittest.TextTestRunner(verbosity=2).run(suite)
