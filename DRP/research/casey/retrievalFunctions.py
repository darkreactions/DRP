import os, sys
full_path = os.path.dirname(os.path.realpath(__file__))+"/"
django_path = full_path[:full_path.rfind("/DRP/")]
if django_path not in sys.path:
  sys.path = [django_path] + sys.path
  os.environ['DJANGO_SETTINGS_MODULE'] = 'DRP.settings'


def get_data_from_ref_file(filename):
  from DRP.models import Data

  with open(django_path+"/"+filename) as f:
    contents = f.read().replace("\n"," ")

    # Create a lower-cased set of refs that we can use for the query.
    refs = contents.split(" ")
    ref_set = set([ref.lower() for ref in refs if ref])

  # Note: does not guarantee that all data will be found.
  data = []
  for datum in Data.objects.all():
    if datum.ref.lower() in ref_set:
      ref_set.remove(datum.ref.lower())
      data.append(datum)

  print ref_set

  return data

