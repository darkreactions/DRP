import os, sys
full_path = os.path.dirname(os.path.realpath(__file__))+"/"
django_path = full_path[:full_path.rfind("/DRP/")]
if django_path not in sys.path:
  sys.path = [django_path] + sys.path
  os.environ['DJANGO_SETTINGS_MODULE'] = 'DRP.settings'


filename = "dataset_2.csv"


def write_uniq_csv(matrix, filename):
  def cleanMatrix(matrix):
    def clean(elem):
      if type(elem)==str or type(elem)==unicode:
        return elem.replace(u"\u2019","'")
      else:
        return elem
    return [[clean(elem) for elem in row] for row in matrix]

  import csv

  prefix = "{}/DRP/research/casey/results/CSVs/".format(django_path)

  with open(prefix + filename, "w") as f:
    writer = csv.writer(f)
    for row in cleanMatrix(matrix):
      try:
        writer.writerow(row)
      except:
        print row
  print "'{}' written!".format(filename)


def split_and_write():
  from DRP.model_building.splitters import default_splitter as splitter
  from DRP.postprocessors import default_postprocessor as postprocessor
  from DRP.preprocessors import default_preprocessor as preprocessor

  from DRP.model_building.generate_models import research_data_filter
  from DRP.retrievalFunctions import get_valid_data
  from DRP.model_building.rxn_calculator import headers
  from DRP.models import Recommendation, Data

  import time

  num_splits = 10

  init_data = get_valid_data()

  orig_data = [headers] + list(research_data_filter(init_data))
  data = preprocessor(orig_data)

  print "All: {}".format(Data.objects.count())
  print "Filtered: {}".format(len(research_data_filter(Data.objects.all())))

  print "All Valid: {}".format(len(get_valid_data()))
  print "Filtered Valid: {}".format(len(orig_data))

  print
  print "Valid..."
  print "Te: {}".format(len([row for row in orig_data if "Te" in row.atoms]))
  print "VTe: {}".format(len([row for row in orig_data if "Te" in row.atoms and "V" in row.atoms]))

  print "Se: {}".format(len([row for row in orig_data if "Se" in row.atoms]))
  print "VSe: {}".format(len([row for row in orig_data if "Se" in row.atoms and "V" in row.atoms]))

  print "Neither: {}".format(len([row for row in orig_data if (not ("Se" in row.atoms and "V" in row.atoms)) and (not ("Te" in row.atoms and "V" in row.atoms))]))

  # Add the Se and Te data.
  headers = data.pop(0) + ["VTe?", "VSe?"]
  data = [row + ["Te" in orig_data[i].atoms and "V" in orig_data[i].atoms,
                 "Se" in orig_data[i].atoms and "V" in orig_data[i].atoms]
          for i, row in enumerate(data)]
  data = [headers] + data


  """
  print "Calculating test_supplement"
  nonsense_recs = Recommendation.objects.filter(nonsense=True)
  pre_nonsense_recs = preprocessor([headers]+list(nonsense_recs))
  nonsense_recs, _ = postprocessor({"all":pre_nonsense_recs}, headers)
  # Change nonsense recs to 0-outcome:
  nonsense_recs = [row[:-1] + [0] for row in nonsense_recs["all"]]
  test_supplement = nonsense_recs[1:] # Remove the "headers"

  data = data + pre_nonsense_recs[1:]
  """

  write_uniq_csv(data, filename)

  """
  for i in xrange(num_splits):
    split_data = splitter(data, headers=headers)
    split_data, new_headers = postprocessor(split_data, headers)

    timestamp = int(time.time())
    write_uniq_csv([new_headers]+split_data["train"],
                   "train_{}.csv".format(timestamp))
    write_uniq_csv([new_headers]+split_data["test"]+test_supplement,
                   "test_{}.csv".format(timestamp))
  """

if __name__=="__main__":
  split_and_write()
