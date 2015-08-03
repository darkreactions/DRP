import unittest
import ConfirmationForm.suite
import LabGroupForm.suite
import LicenseAgreementForm.suite

suite = unittest.TestSuite([
      LabGroupForm.suite,
      ConfirmationForm.suite,
      LicenseAgreementForm.suite
      ])
