from django.core.management.base import BaseCommand, CommandError
from DRP.models import refresh_compound_guide

class Command(BaseCommand):
  def handle(self, *args, **kwargs):
   self.stdout.write("Refresh started!")
   refresh_compound_guide(verbose=True)
   self.stdout.write("Refresh complete!")

