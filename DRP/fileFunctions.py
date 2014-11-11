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


def writeExpandedCSV(name):
  def cleanData(matrix):
    def clean(elem):
      if type(elem)==str or type(elem)==unicode:
        return elem.replace(u"\u2019","")
      else:
        return elem
    return [[clean(elem) for elem in row] for row in matrix]

  from DRP.retrievalFunctions import expand_data,get_expanded_headers,get_valid_data

  print "Gathering data..."
  headers = get_expanded_headers() + ["Lab", "Date Entered"]
  data = expand_data(get_valid_data(),include_lab_info=True, make_new=True, debug=True)

  print "Writing data..."
  data = cleanData(data)

  writeCSV(name, data, headers=headers)

  print "Write complete!"


def createDirIfNecessary(directory):
  if not os.path.exists(directory):
    os.makedirs(directory)
