import sys, os
django_dir = os.path.dirname(os.path.realpath(__file__)).split("DRP")[0]
django_path = "{}/DRP".format(django_dir)
if django_path not in sys.path:
  sys.path.append("{}/DRP".format(django_dir))

os.environ['DJANGO_SETTINGS_MODULE'] = 'DRP.settings'


def writeCSV(name, data, headers=[]):
  import csv

  if headers: data = [headers] + data

  with open(name, "wb") as f:
    csvWriter = csv.writer(f)
    csvWriter.writerows(data)

def writeExpandedCSV(name):
  from DRP.retrievalFunctions import expand_data,get_expanded_headers,get_valid_data
  print "Gathering data..."
  headers = get_expanded_headers() + ["Lab", "Date Entered"]
  data = expand_data(get_valid_data(),include_lab_info=True)

  print "Writing data..."
  writeCSV(name, data, headers=headers)

  print "Write complete!"

