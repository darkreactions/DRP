"""Encapsulates unit tests for various forms."""

import unittest
import CompoundAdminForm
import CompoundCreationForm
import CompoundEditingForm
import ConfirmationForm
import LabGroupForm
import LabGroupJoiningForm
import LicenseAgreementForm

suite = unittest.TestSuite([
    CompoundAdminForm.suite,
    CompoundCreationForm.suite,
    CompoundEditingForm.suite,
    ConfirmationForm.suite,
    LabGroupForm.suite,
    LabGroupJoiningForm.suite,
    LicenseAgreementForm.suite,
])
