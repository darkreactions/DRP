import unittest
import AdvancedFilterForms
import CompoundAdminForm
import CompoundCreationForm
import CompoundEditingForm
import CompoundFilterForm
import CompoundFilterFormSet
import ConfirmationForm
import LabGroupForm
import LabGroupJoiningForm
import LicenseAgreementForm

suite = unittest.TestSuite([
      CompoundAdminForm.suite,
      CompoundCreationForm.suite,
      CompoundEditingForm.suite,
      CompoundFilterForm.suite,
      CompoundFilterFormSet.suite,
      ConfirmationForm.suite,
      LabGroupForm.suite,
      LabGroupJoiningForm.suite,
      LicenseAgreementForm.suite,
      AdvancedFilterForms.suite
      ])
