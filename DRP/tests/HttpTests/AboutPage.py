#!/usr/bin/env python
"""This module contains tests for the about page."""

import unittest
import requests
from HttpTest import GetHttpTest
loadTests = unittest.TestLoader().loadTestsFromTestCase


class AboutPage(GetHttpTest):

    """Performs a get request on the AboutPage to check for Html Validity."""

suite = unittest.TestSuite([
    loadTests(AboutPage)
])

if __name__ == '__main__':
    # Runs the test- a good way to check that this particular test set works without having to run all the tests.
    unittest.main()
