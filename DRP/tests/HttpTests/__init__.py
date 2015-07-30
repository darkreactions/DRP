'''This package provides tests for the responses from pages and views in the DRP project'''
import HomePage.suite
import ContactPage.suite
import AboutPage.suite
import LoginPage.suite
import RegisterPage.suite
import unittest

suite = unittest.TestSuite({
  HomePage.suite,
  ContactPage.suite,
  AboutPage.suite,
  LoginPage.suite,
  RegisterPage.suite
])
