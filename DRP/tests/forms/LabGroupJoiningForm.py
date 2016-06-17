#!/usr/bin/env python
"""The unit test for the LabGroupForm class.
These tests assume that presence tests for teh form fields work as expected
This test suite assumes that Django's check_password functions as
expected.
This test suite assumes that Django's save methods on forms work for
trivial cases.
"""

import unittest
from DRP.tests import DRPTestCase
from DRP.forms import LabGroupJoiningForm
from DRP.models import LabGroup
from django.conf import settings
from django.contrib.auth.models import User
loadTests = unittest.TestLoader().loadTestsFromTestCase


class LegacyPassword(DRPTestCase):

    """Tests for a legacy password check set but no new password has been entered."""

    def setUp(self):
        """Tests for a legacy password check set but no new password has been entered."""
        self.labGroup = LabGroup(title="LegacyPassTest1", address='1, war drobe, Narnia', email='aslan@example.com', legacy_access_code='old_magic')
        self.labGroup.save()
        formData = {'labGroup': self.labGroup.id, 'accessCode': 'old_magic'}
        self.form = LabGroupJoiningForm(formData)

    def test_validation(self):
        """Validates test."""
        valid = self.form.is_valid()
        errString = ''
        for e, m in self.form.errors.items():
            errString += '{0}: {1}\n'.format(e, m)
        self.assertTrue(valid, errString)

    def tearDown(self):
        """Removes lab group created for this test."""
        self.labGroup.delete()


class Password(DRPTestCase):

    """Test that legacy passwords work for joining lab groups."""

    def setUp(self):
        """Create user and labgroupjoiningform for a new user."""
        self.labGroup = LabGroup.objects.makeLabGroup(title="LegacyPassTest1", address='1, war drobe, Narnia', email='aslan@example.com', access_code='old_magic')
        self.labGroup.save()
        formData = {'labGroup': self.labGroup.id, 'accessCode': 'old_magic'}
        self.form = LabGroupJoiningForm(formData)

    def test_validation(self):
        """Validate that the user can join using the form."""
        valid = self.form.is_valid()
        errString = ''
        for e, m in self.form.errors.items():
            errString += '{0}: {1}\n'.format(e, m)
        self.assertTrue(valid, errString)

    def tearDown(self):
        """Delete test lab group."""
        self.labGroup.delete()

suite = unittest.TestSuite([
    loadTests(LegacyPassword),
    loadTests(Password)
])

if __name__ == '__main__':
    unittest.TextTestRunner(verbosity=2).run(suite)
