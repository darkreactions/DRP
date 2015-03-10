import os, sys
full_path = os.path.dirname(os.path.realpath(__file__))+"/"
django_path = full_path[:full_path.rfind("/DRP/")]
if django_path not in sys.path:
  sys.path = [django_path] + sys.path
  os.environ['DJANGO_SETTINGS_MODULE'] = 'DRP.settings'


def get_data_from_ref_file(filename):
  from DRP.retrievalFunctions import get_valid_data

  with open(django_path+"/"+filename) as f:
    contents = f.read()
    refs = contents.split(" ")

    ref_set = set(refs)

  data = get_valid_data()

  data = filter(lambda datum: datum.ref in ref_set, data)

  return data

