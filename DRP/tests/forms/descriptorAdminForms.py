"""Test for the clean method on compound administration forms."""
import unittest
from .baseFormTest import BaseFormTest
from django.conf import settings
from django.contrib.auth.models import User
from DRP.tests.decorators import createsCatRxnDescriptor
from DRP.models import CatRxnDescriptor
from DRP.forms import CatRxnDescriptorForm, CatDescPermittedValueForm

loadTests = unittest.TestLoader().loadTestsFromTestCase


@createsCatRxnDescriptor('family')
class CatCleanTest(BaseFormTest):
    """Test for the cleaning methods."""

    def setUp(self):
        """Set up the form data."""
        desc = CatRxnDescriptor.objects.get(heading='family')
        self.formData = {'id': desc.id,
                         'name': 'edited',
                         'heading': desc.heading
                         }
        self.form = CatRxnDescriptorForm(self.formData, instance=desc)

    def test_validation(self):
        """Test that the validation succeeds."""
        self.validationSucceeds()


@createsCatRxnDescriptor('family', 'nonmanual')
class CatCleanTestFail(BaseFormTest):
    """Test for the cleaning methods."""

    def setUp(self):
        """Set up the form data."""
        desc = CatRxnDescriptor.objects.get(heading='family')
        self.formData = {
            'name': 'edited',
                    'heading': desc.heading
        }
        self.form = CatRxnDescriptorForm(self.formData, instance=desc)


@createsCatRxnDescriptor('family')
class CatValCleanTest(BaseFormTest):
    """Test for the cleaning methods."""

    def setUp(self):
        """Set up the form data."""
        self.formData = {'descriptor': CatRxnDescriptor.objects.get(heading='family').id,
                         'value': 'added'
                         }
        self.form = CatDescPermittedValueForm(self.formData)

    def test_validation(self):
        """Test that the validation succeeds."""
        self.validationSucceeds()


@createsCatRxnDescriptor('family', 'nonmanual')
class CatValCleanTestFail(BaseFormTest):
    """Test for the cleaning methods."""

    def setUp(self):
        """Set up the form data."""
        self.formData = {'descriptor': CatRxnDescriptor.objects.get(heading='family').id,
                         'value': 'added'
                         }
        self.form = CatDescPermittedValueForm(self.formData)

suite = unittest.TestSuite([
    loadTests(CatCleanTest),
    loadTests(CatCleanTestFail),
    loadTests(CatValCleanTest),
    loadTests(CatValCleanTestFail)
])

if __name__ == '__main__':
    unittest.TextTestRunner(verbosity=2).run(suite)
