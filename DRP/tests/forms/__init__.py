"""Encapsulates unit tests for various forms."""

import unittest
from . import compoundAdminForm
from . import compoundCreationForm
from . import compoundEditingForm
from . import confirmationForm
from . import labGroupForm
from . import labGroupJoiningForm
from . import licenseAgreementForm
from . import descriptorAdminForms

suite = unittest.TestSuite([
    compoundAdminForm.suite,
    compoundCreationForm.suite,
    compoundEditingForm.suite,
    confirmationForm.suite,
    labGroupForm.suite,
    labGroupJoiningForm.suite,
    licenseAgreementForm.suite,
    descriptorAdminForms.suite,
])
