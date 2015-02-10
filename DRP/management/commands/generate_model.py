from django.core.management.base import BaseCommand
from DRP.model_building.model_methods import gen_model

class Command(BaseCommand):
  option_list = BaseCommand.option_list

  def handle(self, *args, **options):
    valid_options = ["inactive", "debug", "library"]

    if 2<= len(args) <= 2+len(valid_options):
      #Translate arguments
      self.stdout.write("-- Generating model!...")

      active = not "inactive" in args
      debug = "debug" in args
      library = "sklearn" if "sklearn" in args else "weka"
      model_type = "random_forest"

      title = args[0]
      description = args[1]

      gen_model(title, description, active=active, debug=debug, library=library,
                                    model_type=model_type)

      self.stdout.write("-- Model generation complete.")
    else:
      self.stdout.write("\n--Oops! Check your syntax.")
      self.stdout.write("--Usage: python manage.py generate_model \"Model Name\" \"A brief description of the model.\" [inactive] [debug] [sklearn] [model-type]")

