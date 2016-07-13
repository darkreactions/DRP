"""Test for the clean method on compound administration forms."""
import unittest
from .baseFormTest import BaseFormTest
from django.conf import settings
from django.contrib.auth.models import User
from DRP.tests.decorators import 
from DRP.models import CatRxnDescriptor
from DRP.forms import CatRxnDescriptorForm, CatDescPermittedValueForm 

loadTests = unittest.TestLoader().loadTestsFromTestCase

@createsCatRxnDescriptor('family')
class CatCleanTest(BaseFormTest):
    """Test for the cleaning methods."""

    def setUp(self):
        """Set up the form data."""
        formData = {'id': CatRxnDescriptor.objects.get(heading='family').id
                            'name': 'edited'
                        }
        self.form = CatRxnDescriptorForm(formData) 

    def test_validation(self):
        """Test that the validation succeeds."""
        self.validationSucceeds()

@createsCatRxnDescriptor('family', 'nonmanual')
class CatCleanTestFail(BaseFormTest):
    """Test for the cleaning methods."""

    def setUp(self):
        """Set up the form data."""
        formData = {'id': CatRxnDescriptor.objects.get(heading='family').id,
                            'name': 'edited'
                        }
        self.form = CatRxnDescriptorForm(formData) 


@createsCatRxnDescriptor('family')
class CatValCleanTest(BaseFormTest):
    """Test for the cleaning methods."""

    def setUp(self):
        """Set up the form data."""
        formData = {'descriptor': CatRxnDescriptor.objects.get(heading='family').id
                            'value': 'added'
                        }
        self.form = CatRxnDescPermittedValueForm(formData) 

    def test_validation(self):
        """Test that the validation succeeds."""
        self.validationSucceeds()

@createsCatRxnDescriptor('family', 'nonmanual')
class CatValCleanTestFail(BaseFormTest):
    """Test for the cleaning methods."""

    def setUp(self):
        formData = {'descriptor': CatRxnDescriptor.objects.get(heading='family').id
                            'value': 'added'
                        }
        self.form = CatRxnDescPermittedValueForm(formData) 

suite = unittest.TestSuite([
    loadTests(CatCleanTest),
    loadTests(CatCleanTestFail),
    loadTests(CatValCleanTest),
    loadTests(CatValCleanTestFail)
])

if __name__ == '__main__':
    unittest.TextTestRunner(verbosity=2).run(suite)
