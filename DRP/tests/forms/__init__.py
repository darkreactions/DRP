import unittest
import CompoundAdminForm
import CompoundCreationForm
import CompoundEditingForm
import CompoundFilterForm
import ConfirmationForm
import LabGroupForm
import LabGroupJoiningForm
import LicenseAgreementForm

suite = unittest.TestSuite([
      CompoundAdminForm.suite,
      CompoundCreationForm.suite,
      CompoundEditingForm.suite,
      CompoundFilterForm.suite,
      ConfirmationForm.suite,
      LabGroupForm.suite,
      LabGroupJoiningForm.suite,
      LicenseAgreementForm.suite
      ])
