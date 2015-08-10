#!/usr/bin/env python
'''This module contains tests for teh confirmation page'''

from HttpTests import GetHttpTest, PostHttpTest, OneRedirectionMixin, logsInAs, signsExampleLicense, usesCsrf
from django.contrib.auth.models import User
from DRP.models import ConfirmationCode
from uuid import uuid4
import requests
import unittest

loadTests = unittest.TestLoader().loadTestsFromTestCase

newCompoundUrl = GetHttpRequest.baseUrl + reverse('newCompound')

@logsInAs('Aslan', 'old_magic')
class LicenseRedirect(GetHttpRequest, OneRedirectionMixin):
  '''Tests that the request is redirected if a user tries to view the compound add page in without having
  signed an EULA.'''

  url=newCompoundUrl
  testCodes = ['c9e46ba1-cd2a-4080-88b5-97415fa7c484']

@logsInAs('Aslan', 'old_magic')
@signsExampleLicense('Aslan')
class Lab403Test(GetHttpRequest):
  '''Tests that the view returns the special 403 page when trying to look at a
  compound guide without being in a research group
  '''

  url=newCompoundUrl
  status = 403
  testCodes=['91b3d85b-f975-45f1-b0b5-455d475cfa30']

 

suite = unittest.TestSuite([
  loadTests(ConfirmationPage),
  loadTests(ConfirmationPage2),
  loadTests(PostConfirmationPage),
  loadTests(PostConfirmationPage2)
])

if __name__ == '__main__':
  unittest.main()
