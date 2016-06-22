"""This package provides tests for the responses from pages and views in the DRP project."""
import AboutPage
import AccountPage
import ConfirmationPage
import ContactPage
import HomePage
import LabGroupPage
import LicensePage
import LoginPage
import RegisterPage
import unittest
import CompoundPages
import ReactionPages

suite = unittest.TestSuite([
    AboutPage.suite,
    AccountPage.suite,
    ConfirmationPage.suite,
    ContactPage.suite,
    CompoundPages.suite,
    HomePage.suite,
    LabGroupPage.suite,
    LicensePage.suite,
    LoginPage.suite,
    RegisterPage.suite,
    ReactionPages.suite
])
