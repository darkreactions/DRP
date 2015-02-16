from django.core.management.base import BaseCommand
from DRP.model_building.generate_models import learning_curve
import time

class Command(BaseCommand):
  option_list = BaseCommand.option_list

  def handle(self, *args, **options):
    if 3 <= len(args):
      debug = "debug" in args

      name = args[0]
      description = args[1]
      tags = args[2]

      if "datetime" in args:
        title += "_{}".format(int(time.time()))

      learning_curve(name, description, tags, debug=debug)

    else:
      self.stdout.write("\n--Oops! Check your syntax.")
      self.stdout.write("--Usage: python manage.py learning_curve  \"Model Name Prefix\" \"A brief description of the pipeline.\" \"Tags\" [debug]")

