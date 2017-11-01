"Run descriptor calculations"
import os
import django
import djangofrom DRP import settings
from django.core.management.base import BaseCommand
from DRP.models import Compound, Reaction, NumMolDescriptorValue
from django import db
from django.conf import settings
import logging
import importlib
from django.db import transaction

os.environ['DJANGO_SETTINGS_MODULE'] = 'DRP.settings'
django.setup()

print(Compound.objects.all())

molDescriptorPlugins = [importlib.import_module(plugin) for plugin in settings.MOL_DESCRIPTOR_PLUGINS]

rxnDescriptorPlugins = [importlib.import_module(plugin) for plugin in settings.RXN_DESCRIPTOR_PLUGINS]

print('descriptor plugs \n', molDescriptorPlugins)
print('compounds \n', Compound.objects.all(), '\n\n\n')

for comp in Compound.objects.all():
    print(comp)
    for plugin in molDescriptorPlugins:
        plugin.calculate(comp)
