'''A command to recalculate all descriptors for all compounds when a new plugin has been installed'''
from django.core.management.base import BaseCommand, CommandError
from django.conf import settings
from DRP.models import Compound

descriptorPlugins = [importlib.import_module(plugin) for plugin in settings.MOL_DESCRIPTOR_PLUGINS]

class Command(BaseCommand):
  
  def handle(self, *args, **options):
    for plugin in descriptorPlugins:
      for compound in Compound.objects.all()
        plugin.calculate(compound)
