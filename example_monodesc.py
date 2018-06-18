"""Recalculate the descriptors for all compounds and reactions."""
import os
os.environ['DJANGO_SETTINGS_MODULE'] = 'DRP.settings'
import django

django.setup()




from DRP import settings
from django.core.management.base import BaseCommand
from DRP.models import Compound, Reaction, NumMolDescriptorValue, CompoundQuantity
from django import db
from django.conf import settings
import logging
import importlib
from django.db import transaction

molDescriptorPlugins = [importlib.import_module(plugin) for
			  plugin in settings.MOL_DESCRIPTOR_PLUGINS] # if 'drp' not in str(plugin)]


print('descriptor plugs \n', molDescriptorPlugins)
print('compounds \n', Compound.objects.all(), '\n\n\n')

#NumMolDescriptorValue.objects.all().delete()

# compounds = set([compq.compound for compq in CompoundQuantity.objects.all()])
# compounds = [comp for comp in list(compounds) if NumMolDescriptorValue.objects.filter(compound=comp).count()




counter = 0
for comp in Compound.objects.all():
   counter += 1
   print(comp.smiles)
   print(counter)
   for plugin in molDescriptorPlugins:
       print(NumMolDescriptorValue.objects.filter(compound=comp).count())
       plugin.calculate(comp)
       print(NumMolDescriptorValue.objects.filter(compound=comp).count())
