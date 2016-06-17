#!/usr/bin/env python
"""A module containing tests for the Confirmation Form class.
Since the only novel component is the check for active
and inactive users, this is the only code which is tested.
"""

import unittest
from DRP.tests import DRPTestCase
from DRP.forms import ConfirmationForm
from DRP.models import ConfirmationCode
from django.contrib.auth.models import User

loadTests = unittest.TestLoader().loadTestsFromTestCase

class InactiveUser(DRPTestCase):
    """Tests the case where good credentials for an inactive user have been submitted."""

    def setUp(self):
        """Set up."""
        self.username = 'Aslan'
        self.password = 'old_magic'
        self.user = User.objects.create_user(username=self.username, password=self.password)
        self.user.is_active = False
        self.user.save()

    def runTest(self):
        """Run test."""
        form = ConfirmationForm(data={'username': self.username, 'password': self.password})
        valid = form.is_valid()
        errString = ''
        for e, m in form.errors.items():
            errString += '{0}: {1}\n'.format(e, m)
        self.assertTrue(valid, errString)

    def tearDown(self):
        """Delete user object created for the test."""
        self.user.delete()


class ActiveUser(DRPTestCase):

    """Tests the case where good credentials for an active user have been sumbitted."""

    def setUp(self):
        """Set up a user for the test."""
        self.username = 'Aslan'
        self.password = 'old_magic'
        self.user = User.objects.create_user(username=self.username, password=self.password)
        self.user.is_active = True
        self.user.save()

    def runTest(self):
        """Run test to make sure that the form does not evaluate this as valid."""
        form = ConfirmationForm(data={'username': self.username, 'password': self.password})
        self.assertFalse(form.is_valid())

    def tearDown(self):
        """Delete user object created for the test."""
        self.user.delete()

suite = unittest.TestSuite([
    loadTests(InactiveUser),
    loadTests(ActiveUser)
])

if __name__ == '__main__':
    # Runs the test- a good way to check that this particular test set works without having to run all the tests.
    unittest.main()
