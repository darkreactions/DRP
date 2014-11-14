from django.core.management.base import BaseCommand, CommandError
from DRP.models import refresh_compound_guide

class Command(BaseCommand):
  def handle(self, *args, **kwargs):
   self.stdout.write("Refresh started!")

   debug = "debug" in args
   clear = "clear" in args
   refresh_compound_guide(verbose=True, debug=debug, clear=clear)

   self.stdout.write("Refresh complete!")

