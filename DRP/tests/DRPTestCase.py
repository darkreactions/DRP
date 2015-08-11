'''A module containing only the DRPTestCase class'''
import unittest
from django.conf import settings
import DRP.models as models 
from django.contrib.auth.models import User

class DRPTestCase(unittest.TestCase):
  '''A quick and dirty safety valve to stop people accidentally running database tests in production environments
  For more information see the documentation.
  '''

  def __init__(self, *args, **kwargs):
    if not settings.TESTING:
      raise RuntimeError('Testing environment not enabled')
    else:
      super(DRPTestCase, self).__init__(*args, **kwargs)
      self.addCleanup(cleanUpDatabase)

def cleanUpDatabase():
  models.Compound.objects.all().delete()
  models.ChemicalClass.objects.all().delete()
  models.ConfirmationCode.objects.all().delete()
  models.DataSet.objects.all().delete()
  models.LabGroup.objects.all().delete()
  models.LicenseAgreement.objects.all()
  models.License.objects.all().delete()
  models.LegacyStatsModel.objects.all().delete()
  models.StatsModel.objects.all().delete()
  models.PerformedReaction.objects.all().delete()
  models.MolDescriptor.objects.all().delete()
  models.RxnDescriptor.objects.all().delete()
  models.MolDescriptorValue.objects.all().delete()
  models.RxnDescriptorValue.objects.all().delete()
  models.StatsModelTag.objects.all().delete()
  User.objects.all().exclude(username='root').delete()

def runTests(suite):
  '''A function which empties out the database prior to and after running the tests contained in suite'''
  cleanUpDatabase()
  unittest.TextTestRunner(verbosity=2).run(suite)
  cleanUpDatabase()
