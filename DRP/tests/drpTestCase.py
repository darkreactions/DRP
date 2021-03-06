"""A module containing only the DRPTestCase class."""
import unittest
from django.conf import settings
import importlib
import DRP
from django.contrib.auth.models import User
import logging
tracer = logging.getLogger('DRP.tracer')
# this prevents a cyclic dependency problem
molDescriptorPlugins = [importlib.import_module(
    plugin) for plugin in settings.MOL_DESCRIPTOR_PLUGINS]
# this prevents a cyclic dependency problem
rxnDescriptorPlugins = [importlib.import_module(
    plugin) for plugin in settings.RXN_DESCRIPTOR_PLUGINS]


class DRPTestCase(unittest.TestCase):
    """
    A quick and dirty safety valve to stop people accidentally running database tests in production environments.

    For more information see the documentation.
    """

    deleteDescriptors = True

    def __init__(self, *args, **kwargs):
        """Break start of a test if testing mode is not set."""
        self.cleaned = False
        """Quick break if it looks like you're about to really delete things."""
        if not settings.TESTING:
            raise RuntimeError('Testing environment not enabled')
        else:
            super(DRPTestCase, self).__init__(*args, **kwargs)

    def __getattribute__(self, attribute_name):
        """Clear the database if setting up a test."""
        if attribute_name == 'setUp' and self.cleaned is False:
            cleanUpDatabase(self.deleteDescriptors)
            self.cleaned = True
        return super(DRPTestCase, self).__getattribute__(attribute_name)

    def doCleanups(self):
        """Clean the database if required."""
        self.cleaned = False
        return super(DRPTestCase, self).doCleanups()


def cleanUpDatabase(deleteDescriptors=False):
    """Truncate the whole database."""
    tracer.debug('Erasing DB')
    DRP.models.CompoundQuantity.objects.all().delete()
    DRP.models.Compound.objects.all().delete()
    DRP.models.CompoundRole.objects.all().delete()
    DRP.models.ChemicalClass.objects.all().delete()
    DRP.models.ConfirmationCode.objects.all().delete()
    DRP.models.LicenseAgreement.objects.all()
    DRP.models.License.objects.all().delete()
    DRP.models.StatsModel.objects.all().delete()
    DRP.models.DataSetRelation.objects.all().delete()
    DRP.models.PerformedReaction.objects.all().delete()
    DRP.models.CatMolDescriptorValue.objects.all().delete()
    DRP.models.BoolMolDescriptorValue.objects.all().delete()
    DRP.models.NumMolDescriptorValue.objects.all().delete()
    DRP.models.OrdMolDescriptorValue.objects.all().delete()
    DRP.models.LabGroup.objects.all().delete()
    if deleteDescriptors:
        keepDescriptors = []
        for plugin in molDescriptorPlugins:
            keepDescriptors += [plugin.descriptorDict[key]
                                .descriptor_ptr.pk for key in plugin.descriptorDict]
        for plugin in rxnDescriptorPlugins:
            keepDescriptors += [plugin.descriptorDict[key]
                                .descriptor_ptr.pk for key in plugin.descriptorDict]
        DRP.models.Descriptor.objects.exclude(pk__in=keepDescriptors).delete()
        User.objects.all().exclude(username='root').delete()


def runTests(suite, failfast=False):
    """A function which empties out the database prior to and after running the tests contained in suite."""
    return unittest.TextTestRunner(verbosity=4, failfast=failfast).run(suite)
