#!/usr/bin/env python

import unittest
import os

from DRP.models import BoolRxnDescriptor, PerformedReaction, ModelContainer, Descriptor, NumericDescriptor
from django.conf import settings
loadTests = unittest.TestLoader().loadTestsFromTestCase
from DRPTestCase import DRPTestCase, runTests
from DRP.recommender import SampleGridParameters, ReactionSieve, ReactionRecommender




class ModelTest(DRPTestCase):
    """Tests basic recommender functions."""

    def setup(self):
        model_container = ModelContainer.objects.get(pk=1) # TODO add more general model


        desc_dict = {}
        desc_dict[NumericDescriptor.objects.get(pk=6)] = [1, 2, 3, 8]
        desc_dict[NumericDescriptor.objects.get(pk=5)] = [1440, 3600]
        desc_dict[NumericDescriptor.objects.get(pk=4)] = [363.15, 423.15]

        desired_desc_dict = {BoolRxnDescriptor.objects.get(pk=2):[True]}

        self.recommender = None


suite = unittest.TestSuite([loadTests(ModelTest)])

if __name__=='__main__':
    os.environ['DJANGO_SETTINGS_MODULE'] = 'DRP.settings'
    import django
    django.setup()
    runTests(suite)
    print 'okay'