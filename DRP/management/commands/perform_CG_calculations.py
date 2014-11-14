from django.core.management.base import BaseCommand, CommandError
from optparse import make_option
from DRP.models import perform_CG_calculations

class Command(BaseCommand):
  option_list = BaseCommand.option_list + (
  make_option('--only-missing',
       action='store_true',
       dest='only_missing',
       default=False,
       help='Perform calculations on only the missing data.'),
  ) + (
  make_option('--attempt-failed',
       action='store_true',
       dest='attempt_failed',
       default=True,
       help='Perform calculations on previously failed calculations.'),
  )

  def handle(self, *args, **options):
   #Translate arguments
   self.stdout.write("Calculations started!")
   perform_CG_calculations(only_missing=options["only_missing"], attempt_failed=options["attempt_failed"],verbose=True)
   self.stdout.write("CG_calculations construction complete!")

