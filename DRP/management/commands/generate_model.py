from django.core.management.base import BaseCommand, CommandError
from DRP.model_building.model_methods import gen_model

class Command(BaseCommand):
  option_list = BaseCommand.option_list

  def handle(self, *args, **options):
    if len(args)==2:
      #Translate arguments
      self.stdout.write("-- Generating model!...")
      title = args[0]
      description = args[1]
      gen_model(title, description)
      self.stdout.write("-- Model generation complete.")
    else:
      self.stdout.write("\n--Oops! Check your syntax.")
      self.stdout.write("--Usage: python manage.py generate_model \"Model Name\" \"A brief description of the model.\"")

