from django.core.management.base import BaseCommand, CommandError
from DRP.update import update_data_calcs

class Command(BaseCommand):
  option_list = BaseCommand.option_list

  def handle(self, *args, **options):
    if len(args)==1:

      debug = True

      update_data_calcs()

      self.stdout.write("-- Calculations updated!")
    else:
      self.stdout.write("--Oops! Can't understand input...")
      self.stdout.write("--Usage: python manage.py update_data_calcs)
