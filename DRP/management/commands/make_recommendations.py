from django.core.management.base import BaseCommand, CommandError
from DRP.recommendation.recommend import create_new_recommendations
from DRP.models import get_Lab_Group

class Command(BaseCommand):
  option_list = BaseCommand.option_list

  def handle(self, *args, **options):
    if len(args)==1:

      debug = True

      lab_group = get_Lab_Group(args[0])
      create_new_recommendations(lab_group, debug=debug, bare_debug=debug)

      self.stdout.write("-- Recommendations complete!")
    else:
      self.stdout.write("--Oops! Can't understand input...")
      self.stdout.write("--Usage: python manage.py make_recommendations \"Lab Title\"")
