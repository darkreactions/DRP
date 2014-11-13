from django.core.management.base import BaseCommand, CommandError
from DRP.update import update_data_calcs

class Command(BaseCommand):
  def handle(self, *args, **options):
    if len(args)<=2:

      debug = "debug" in args
      clear = "clear" in args
      update_data_calcs(debug=debug, clear=clear)

      self.stdout.write("-- Calculations updated!")
    else:
      self.stdout.write("--Oops! Can't understand input...")
      self.stdout.write("--Usage: python manage.py update_data_calcs [debug]")
