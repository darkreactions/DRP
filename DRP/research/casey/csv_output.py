import os, sys
full_path = os.path.dirname(os.path.realpath(__file__))+"/"
django_path = full_path[:full_path.rfind("/DRP/")]
if django_path not in sys.path:
  sys.path = [django_path] + sys.path
  os.environ['DJANGO_SETTINGS_MODULE'] = 'DRP.settings'

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
  from DRP.models import Recommendation

  from django.db.models import Q
  import time

  num_splits = 10

  data = get_valid_data()

  data = [headers] + list(research_data_filter(data))
  data = preprocessor(data)

  print "Calculating train_supplement"
  calculated = ~Q(calculations=None)
  nonsense_recs = Recommendation.objects.filter(Q(nonsense=True), calculated)

  clean_recs = []
  for i, rec in enumerate(nonsense_recs):

    try:
      rec.get_calculations_list()
      clean_recs.append(rec)
    except Exception as e:
      pass

  nonsense_recs = preprocessor([headers]+clean_recs)
  nonsense_recs, _ = postprocessor({"all":nonsense_recs}, headers)
  train_supplement = nonsense_recs["all"][1:] # Ignore the duplicate headers
  print len(train_supplement)

  write_uniq_csv(data, "all.csv")

  for i in xrange(num_splits):
    split_data = splitter(data, headers=headers)
    split_data, new_headers = postprocessor(split_data, headers)

    timestamp = int(time.time())
    write_uniq_csv([new_headers]+split_data["train"]+train_supplement,
                   "train_{}.csv".format(timestamp))
    write_uniq_csv([new_headers]+split_data["test"]+train_supplement,
                   "test_{}.csv".format(timestamp))


if __name__=="__main__":
  split_and_write()
