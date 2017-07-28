"""This package provides tests for the responses from pages and views in the DRP project."""
from . import aboutPage
from . import accountPage
from . import confirmationPage
from . import contactPage
from . import homePage
from . import labGroupPage
from . import licensePage
from . import loginPage
from . import registerPage
import unittest
from . import compoundPages
from . import reactionPages

suite = unittest.TestSuite([
    aboutPage.suite,
    accountPage.suite,
    confirmationPage.suite,
    contactPage.suite,
    compoundPages.suite,
    homePage.suite,
    labGroupPage.suite,
    licensePage.suite,
    loginPage.suite,
    registerPage.suite,
    reactionPages.suite
])
