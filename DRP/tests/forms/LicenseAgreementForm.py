#!/usr/bin/env python
'''The tests for the LicenseAgreementForm.
These tests assume that html validation will be covered as a compenent of template testing
These tests assume that the user authentication provided by the Django Form work as expected.
These tests assume that the very simple "save" method works as expected (and in any case will be tested by the view tests)
'''

import unittest
from DRP.tests import DRPTestCase
from django.contrib.auth.models import User
from DRP.models import License
from DRP.forms import LicenseAgreementForm
import datetime

loadTests = unittest.TestLoader().loadTestsFromTestCase

class LicenceAgreementForm(DRPTestCase):

  def setUp(self):
    self.passes = ['banana', 'turkishdelight']
    self.user = User.objects.create_user(username='Aslan', password=self.passes[0])
    self.user.save()
    self.license = License(text='This is some text', effectiveDate = datetime.date.today() - datetime.timedelta(1)) #license was created yesterday, always.
    self.license.save()
    self.invUser = User.objects.create_user(username='Whitewitch', password=self.passes[1]) 
    self.invUser.save()
    self.invlicense = License(text='This is an invalid license', effectiveDate = datetime.date.today() + datetime.timedelta(1)) #license was created tomorrow, always.
    self.invlicense.save()
    self.combinations = (
      (self.user, self.passes[0], self.license.id, True), #tuple of the form, user, submitted license id, expected is_valid form value assumes 'current' user is Aslan in all cases.
      (self.user, self.passes[0], self.invlicense.id, False),
      (self.invUser, self.passes[1], self.license, False),
      (self.invUser, self.passes[1], self.invlicense, False)
    )

  def runTest(self):
    i=0
    for c in self.combinations:
      i += 1
      form = LicenseAgreementForm(self.user, self.license, data={'username':c[0].username, 'password':c[1], 'licenseId':c[2]})
      errString = 'User: {0}, Pass:{1}, license_id:{2} (valid={3}); row {4} [{5}, {6}]\n'.format(c[0].username, c[1], c[2], self.license.id, i, form.is_valid(), c[3])
      for k, v in form.errors.items():
        errString += '\n'.join(k + ':' + e for e in v)
      self.assertEqual(c[3], form.is_valid(), errString)
  
  def tearDown(self):
    self.user.delete()
    self.invUser.delete()
    self.license.delete()
    self.invlicense.delete()

suite =  unittest.TestSuite([
  loadTests(LicenseAgreementForm)  
])

if __name__=='__main__':
  unittest.main()
