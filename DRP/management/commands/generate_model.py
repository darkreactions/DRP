from django.core.management.base import BaseCommand, CommandError
from DRP.model_building.model_methods import gen_model

class Command(BaseCommand):
  option_list = BaseCommand.option_list

  def handle(self, *args, **options):
    if 2<= len(args) <=4:
      #Translate arguments
      self.stdout.write("-- Generating model!...")

      active = not "inactive" in args
      debug = "debug" in args
      args = filter(lambda arg: arg!="inactive" and arg!="debug", args)

      title = args[0]
      description = args[1]

      gen_model(title, description, active=active, debug=debug)

      self.stdout.write("-- Model generation complete.")
    else:
      self.stdout.write("\n--Oops! Check your syntax.")
      self.stdout.write("--Usage: python manage.py generate_model \"Model Name\" \"A brief description of the model.\" [inactive] [debug]")

