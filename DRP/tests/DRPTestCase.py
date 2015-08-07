'''A module containing only the DRPTestCase class'''
from unittest import TestCase 
from django.conf import settings
import DRP.models as models 
from django.contrib.auth.models import User

class DRPTestCase(TestCase):
  '''A quick and dirty safety valve to stop people accidentally running database tests in production environments
  For more information see the documentation.
  '''

  def __init__(self, *args, **kwargs):
    if not settings.TESTING:
      raise RuntimeError('Testing environment not enabled')
    else:
      super(DRPTestCase, self).__init__(*args, **kwargs)
      models.Compound.objects.all().delete()
      models.ChemicalClass.objects.all().delete()
      models.LabGroup.objects.all().delete()
      models.ConfirmationCode.objects.all().delete()
      models.LicenseAgreement.objects.all()
      models.License.objects.all().delete()
      User.objects.all().exclude(username='root').delete()
      
      
