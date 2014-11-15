import os, sys
full_path = os.path.dirname(os.path.realpath(__file__))+"/"
django_path = full_path[:full_path.rfind("/DRP/")]
if django_path not in sys.path:
  sys.path = [django_path] + sys.path
  os.environ['DJANGO_SETTINGS_MODULE'] = 'DRP.settings'


def recalculate_usable_models():
  from DRP.models import ModelStats

  # Variable Setup
  all_models = ModelStats.objects.all()
  good = 0

  for model in all_models:
    usable = model.check_usability()
    if usable: good += 1

  print "{} of {} models usable!".format(good, all_models.count())

if __name__=="__main__":
  recalculate_usable_models()

