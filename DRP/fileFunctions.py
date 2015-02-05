import os, sys
full_path = os.path.dirname(os.path.realpath(__file__))+"/"
django_path = full_path[:full_path.rfind("/DRP/")]
if django_path not in sys.path:
  sys.path = [django_path] + sys.path
  os.environ['DJANGO_SETTINGS_MODULE'] = 'DRP.settings'

def writeCSV(name, data, headers=[]):
  import csv

  if headers: data = [headers] + data

  with open(name, "wb") as f:
    csvWriter = csv.writer(f)
    csvWriter.writerows(data)


def writeExpandedCSV(filename, debug=True, include_lab_info=True):
  def cleanData(matrix):
    def clean(elem):
      if type(elem)==str or type(elem)==unicode:
        return elem.replace(u"\u2019","")
      else:
        return elem
    return [[clean(elem) for elem in row] for row in matrix]

  from DRP.retrievalFunctions import expand_data,get_expanded_headers,get_valid_data

  if debug: print "Gathering data..."

  headers = get_expanded_headers()

  if include_lab_info:
    headers.extend(["Lab", "Date Entered"])

  data = expand_data(get_valid_data(), include_lab_info=include_lab_info,
                                       make_new=True,
                                       debug=debug)

  if debug: print "Writing data..."
  data = cleanData(data)

  writeCSV(filename, data, headers=headers)

  if debug: print "Write complete!"


def createDirIfNecessary(directory):
  if not os.path.exists(directory):
    os.makedirs(directory)

def get_django_path():
  import os
  full_path = os.path.dirname(os.path.realpath(__file__))+"/"
  django_path = full_path[:full_path.rfind("/DRP/")]
  return django_path
