"""Encapsulates unit tests for various forms."""

import unittest
import compoundAdminForm
import compoundCreationForm
import compoundEditingForm
import confirmationForm
import labGroupForm
import labGroupJoiningForm
import licenseAgreementForm

suite = unittest.TestSuite([
    compoundAdminForm.suite,
    compoundCreationForm.suite,
    compoundEditingForm.suite,
    confirmationForm.suite,
    labGroupForm.suite,
    labGroupJoiningForm.suite,
    licenseAgreementForm.suite,
])
