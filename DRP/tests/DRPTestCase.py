'''A module containing only the DRPTestCase class'''
from unittest import TestCase 
from django.conf import settings

class DRPTestCase(TestCase):
  '''A quick and dirty safety valve to stop people accidentally running database tests in production environments
  For more information see the documentation.
  '''

  def __init__(self, *args, **kwargs):
    if not settings.TESTING:
      raise RuntimeError('Testing environment not enabled')
    else:
      super(DRPTestCase, self).__init__(*args, **kwargs)
