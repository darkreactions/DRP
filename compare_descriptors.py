## Using the app from inside a script
import os
os.environ['DJANGO_SETTINGS_MODULE'] = 'DRP.settings'
import django

django.setup()

from DRP import settings
from django.core.management.base import BaseCommand
from django import db
from django.conf import settings
import logging
import importlib
from django.db import transaction
from DRP.models import PerformedReaction, Descriptor, NumRxnDescriptorValue, Compound, Reaction, CompoundQuantity


PerformedReaction.objects.all().toCsv(open('test.csv', 'w'), expanded=True) #, labGroups=["Norquist Lab"])
