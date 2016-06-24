"""This package provides tests for the responses from pages and views in the DRP project."""
import .aboutPage
import .accountPage
import .confirmationPage
import .contactPage
import .homePage
import .labGroupPage
import .licensePage
import .loginPage
import .registerPage
import .unittest
import .compoundPages
import .reactionPages

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
