#!/usr/bin/env python
"""Tests for custom model validators."""
import unittest
from .drpTestCase import DRPTestCase, runTests
from DRP.models import validators
from django.core.exceptions import ValidationError
import datetime

loadTests = unittest.TestLoader().loadTestsFromTestCase


class NotInTheFuture(DRPTestCase):
    """Ensures that the validator for things in the future functions properly."""

    def setUp(self):
        """Blank out normal actions."""
        pass

    def test_future(self):
        """Test a date in the future."""
        self.assertRaises(ValidationError, validators.notInTheFuture,
                          datetime.datetime.now() + datetime.timedelta(1))

    def test_past(self):
        """Test a date in the past."""
        validators.notInTheFuture(
            datetime.datetime.now() - datetime.timedelta(1))

    def test_present(self):
        """Test today."""
        validators.notInTheFuture(datetime.datetime.now())

    def tearDown(self):
        """Blank out normal actions."""
        pass


class GreaterThan(DRPTestCase):
    """Tests the greater than decorator."""

    def setUp(self):
        """Set up an instance of the validator."""
        self.validator = validators.GreaterThanValidator(0)

    def test_less(self):
        """Test a less than value."""
        self.assertRaises(ValidationError, self.validator, -1)

    def test_equal(self):
        """Test an equal value."""
        self.assertRaises(ValidationError, self.validator, 0)

    def test_greater(self):
        """Test a greater value."""
        self.validator(2)

    def tearDown(self):
        """Override."""
        pass

suite = unittest.TestSuite([
    loadTests(GreaterThan),
    loadTests(NotInTheFuture)
])

if __name__ == '__main__':
    runTests(suite)
    # Runs the test- a good way to check that this particular test set works
    # without having to run all the tests.
