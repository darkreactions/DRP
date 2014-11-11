from django.core.management.base import BaseCommand, CommandError
from DRP.model_building.model_methods import gen_model

class Command(BaseCommand):
  option_list = BaseCommand.option_list

  def handle(self, *args, **options):
    if len(args)==2 or len(args)==3:
      #Translate arguments
      self.stdout.write("-- Generating model!...")
      active = True
      if len(args)==3:
        active = all(["inactive"!=arg for arg in args])
        args = filter(lambda arg: arg!="inactive", args)

        if len(args)==3:
          raise Exception("Unknown argument used -- perhaps you meant 'inactive'?")

      title = args[0]
      description = args[1]

      gen_model(title, description, active=active)

      self.stdout.write("-- Model generation complete.")
    else:
      self.stdout.write("\n--Oops! Check your syntax.")
      self.stdout.write("--Usage: python manage.py generate_model \"Model Name\" \"A brief description of the model.\" [inactive]")

