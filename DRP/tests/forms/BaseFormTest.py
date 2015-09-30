'''A base form test for testing forms'''

from DRP.tests import DRPTestCase
  
class BaseFormTest(DRPTestCase):

  def setUp(self):
    self.setUpFormData()

  def test_validation(self):
    self.validationFails()

  def validationFails(self):
    '''a test for cases where the form should not validate'''
    self.assertFalse(self.form.is_valid(), 'Form should have failed... submitted data: {0}'.format(self.formData))

  def validationSucceeds(self):
    '''Test that the form does validate'''
    valid = self.form.is_valid()
    errString = ''
    for e, m in self.form.errors.items():
      errString += '{0}: {1}\n'.format(e, m)
    errString += ' Submitted data = {0}'.format(self.formData)
    self.assertTrue(valid, errString)

