'''A module containing only the DRPTestCase class'''
import unittest
from django.conf import settings
import importlib
import DRP
from django.contrib.auth.models import User
molDescriptorPlugins = [importlib.import_module(plugin) for plugin in settings.MOL_DESCRIPTOR_PLUGINS] #this prevents a cyclic dependency problem

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
  DRP.models.CompoundQuantity.objects.all().delete()
  DRP.models.Compound.objects.all().delete()
  DRP.models.CompoundRole.objects.all().delete()
  DRP.models.ChemicalClass.objects.all().delete()
  DRP.models.ConfirmationCode.objects.all().delete()
  DRP.models.LicenseAgreement.objects.all()
  DRP.models.License.objects.all().delete()
  DRP.models.LegacyStatsModel.objects.all().delete()
  DRP.models.StatsModel.objects.all().delete()
#  DRP.models.TestSet.objects.all().delete()
#  DRP.models.TestSetRelation.objects.all().delete()
#  DRP.models.TrainingSet.objects.all().delete()
  DRP.models.PerformedReaction.objects.all().delete()
  DRP.models.CatMolDescriptorValue.objects.all().delete()
  DRP.models.BoolMolDescriptorValue.objects.all().delete()
  DRP.models.NumMolDescriptorValue.objects.all().delete()
  DRP.models.OrdMolDescriptorValue.objects.all().delete()
#  DRP.models.StatsModelTag.objects.all().delete()
  DRP.models.LabGroup.objects.all().delete()
  keepDescriptors = []
  for plugin in molDescriptorPlugins:
    keepDescriptors += [plugin.descriptorDict[key].descriptor_ptr.pk for key in plugin.descriptorDict]
  DRP.models.Descriptor.objects.exclude(pk__in=keepDescriptors).delete()
  User.objects.all().exclude(username='root').delete()

def runTests(suite):
  '''A function which empties out the database prior to and after running the tests contained in suite'''
  cleanUpDatabase()
  return unittest.TextTestRunner(verbosity=4).run(suite)
