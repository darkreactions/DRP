#!/usr/bin/env python

#This file contains a (very) loose framework from which others can base their test files and be
#conformant with the local arrangement of test cases.
#For more information about the structure of tests, consult the python documentation at
#https://docs.python.org/2/library/unittest.html

#REMEMBER when you create a new test case to add it to the suite() method, and then
#to have that suite method called in AllTests.py

import unittest
from DRPTestCase import DRPTestCase, runTests
loadTests = unittest.TestLoader().loadTestsFromTestCase
from DRP.models import validators
from django.core.exceptions import ValidationError
import datetime

class NotInTheFuture(DRPTestCase):
    #This class exemplifies the standard structure of a test. Check the documentation for 'rolling your own'

    def setUp(self):
        pass

    def test_future(self):
        self.assertRaises(ValidationError, validators.notInTheFuture, datetime.datetime.now() + datetime.timedelta(1))

    def test_past(self):
        validators.notInTheFuture(datetime.datetime.now() - datetime.timedelta(1))

    def text_present(self):
        validators.notInTheFuture(datetime.datetime.now())

    def tearDown(self):
        pass

class GreaterThan(DRPTestCase):

    def setUp(self):
        self.validator = validators.GreaterThanValidator(0)

    def test_less(self):
        self.assertRaises(ValidationError, self.validator, -1)

    def test_equal(self):
        self.assertRaises(ValidationError, self.validator, 0)

    def test_greater(self):
        self.validator(2)

    def tearDown(self):
        pass

suite = unittest.TestSuite([
            loadTests(GreaterThan),
            loadTests(NotInTheFuture)
            ])

if __name__=='__main__':
    runTests(suite)
    #Runs the test- a good way to check that this particular test set works without having to run all the tests.
