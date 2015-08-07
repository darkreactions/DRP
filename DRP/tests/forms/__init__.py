import unittest
import CompoundAdminForm
import CompoundCreationForm
import ConfirmationForm
import LabGroupForm
import LabGroupJoiningForm
import LicenseAgreementForm

suite = unittest.TestSuite([
      CompoundAdminForm.suite,
      CompoundCreationForm.suite,
      ConfirmationForm.suite,
      LabGroupForm.suite,
      LabGroupJoiningForm.suite,
      LicenseAgreementForm.suite
      ])
