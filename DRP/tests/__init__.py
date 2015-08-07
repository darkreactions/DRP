from DRPTestCase import DRPTestCase
import unittest
import Email
import forms
import HttpTests

suite = unittest.TestSuite([
  Email.suite,
  forms.suite,
  HttpTests.suit
]) 
