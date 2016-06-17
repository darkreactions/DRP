"""A base form test for testing forms."""

from DRP.tests import DRPTestCase


class BaseFormTest(DRPTestCase):

    """Testing basic form validation."""

    def setUp(self):
        """Set up a form to test validation."""
        self.setUpFormData()

    def test_validation(self):
        """A test for cases where the form should validate."""
        self.validationFails()

    def validationFails(self):
        """A test for cases where the form should not validate."""
        self.assertFalse(self.form.is_valid(), 'Form should have failed... submitted data: {0}'.format(self.formData))

    def validationSucceeds(self):
        """Test that the form does validate."""
        valid = self.form.is_valid()
        errString = ''
        if hasattr(self.form.errors, 'items'):
            for e, m in self.form.errors.items():
                errString += '{0}: {1}\n'.format(e, m)
        else:
            for s in self.form.errors:
                for e, m in s.items():
                    errString += '{0}: {1}\n'.format(e, m)
        errString += ' Submitted data = {0}'.format(self.formData)
        self.assertTrue(valid, errString)
