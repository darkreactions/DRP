from django.core.management.base import BaseCommand, CommandError
from DRP.recommendation.recommend import create_new_recommendations
from DRP.fileFunctions import writeExpandedCSV

class Command(BaseCommand):
  option_list = BaseCommand.option_list

  def handle(self, *args, **options):
    if len(args)==1:

      debug = False
      include_lab_info = False

      writeExpandedCSV(args[0], debug=debug, include_lab_info=include_lab_info)

      self.stdout.write("-- File written!")
    else:
      self.stdout.write("--Oops! Can't understand input...")
      self.stdout.write("--Usage: python manage.py write_expanded_csv\"")

