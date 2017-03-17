"""Recalculate the descriptors for all compounds and reactions."""

import django
django.setup()


from DRP import settings
from django.core.management.base import BaseCommand
from DRP.models import Compound, Reaction, NumMolDescriptorValue
from django import db
from django.conf import settings
import logging
import importlib
from django.db import transaction


molDescriptorPlugins = [importlib.import_module(plugin) for
                                                plugin in settings.MOL_DESCRIPTOR_PLUGINS]

rxnDescriptorPlugins = [importlib.import_module(plugin) for
                                                plugin in settings.RXN_DESCRIPTOR_PLUGINS]

print('descriptor plugs \n', molDescriptorPlugins)
print('compounds \n', Compound.objects.all(), '\n\n\n')

NumMolDescriptorValue.objects.all().delete()


for comp in Compound.objects.all():
    print(comp)
    for plugin in molDescriptorPlugins:
        plugin.calculate(comp)

