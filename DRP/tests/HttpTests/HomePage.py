#!/usr/bin/env python
"""This package provides tests for the home page."""
from HttpTest import GetHttpTest
import unittest
import requests

loadTests = unittest.TestLoader().loadTestsFromTestCase

class HomePage(GetHttpTest):
    """Tests the home page for HTML validity, and that the correct template is rendered."""

    url = GetHttpTest.baseUrl
    testCodes = ['2240b1ff-895c-458f-adf2-d04d85a164d1', 'd68a82db-bd18-4a9f-a1a2-03b3bb259595']
    """A list of unique IDs. A check is performed to see that these are in the returned HTML."""

suite = unittest.TestSuite([
    loadTests(HomePage)
])

if __name__ == '__main__':
    # Runs the test- a good way to check that this particular test set works without having to run all the tests.
    unittest.main()
