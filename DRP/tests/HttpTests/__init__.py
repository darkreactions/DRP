'''This package provides tests for the responses from pages and views in the DRP project'''
import HomePage
import ContactPage
import AboutPage
import LoginPage
import RegisterPage
import LicensePage
import unittest

suite = unittest.TestSuite([
  HomePage.suite,
  ContactPage.suite,
  AboutPage.suite,
  LoginPage.suite,
  RegisterPage.suite
])
