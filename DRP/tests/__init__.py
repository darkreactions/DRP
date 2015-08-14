from DRPTestCase import DRPTestCase, runTests
import unittest
import Email
import forms
import HttpTests
import decorators
import CompoundFromCSV

suite = unittest.TestSuite([
  Email.suite,
  forms.suite,
  HttpTests.suite,
  CompoundFromCSV.suite
]) 
