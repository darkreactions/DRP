from django.core.management.base import BaseCommand
from DRP.model_building.generate_models import gen_model
import time

class Command(BaseCommand):
  option_list = BaseCommand.option_list

  def handle(self, *args, **options):
    valid_options = ["inactive", "debug", "library"]

    if 2<= len(args) <= 2+len(valid_options):
      active = not "inactive" in args
      debug = "debug" in args

      title = args[0]
      description = args[1]

      if "datetime_title" in args:
        title += "_{}".format(int(time.time()))

      gen_model(title, description, active=active, debug=debug)

    else:
      self.stdout.write("\n--Oops! Check your syntax.")
      self.stdout.write("--Usage: python manage.py generate_model \"Model Name\" \"A brief description of the model.\" [inactive] [debug]")

