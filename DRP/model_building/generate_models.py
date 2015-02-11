#!/usr/bin/env python

#Grab the Django settings if they aren't already set.
import os, sys
full_path = os.path.dirname(os.path.realpath(__file__))+"/"
django_path = full_path[:full_path.rfind("/DRP/")]
if django_path not in sys.path:
  sys.path = [django_path] + sys.path
  os.environ['DJANGO_SETTINGS_MODULE'] = 'DRP.settings'

def _default_preprocessor(data):
  def is_num(elem):
    try:
      float(elem)
      return True
    except:
      return False

  def norm(elems):
    elems = map(float,elems)
    minimum = min(elems)
    maximum = max(elems)
    diff = float(maximum-minimum)

    return [(maximum-elem)/diff for elem in elems]

  # Variable Setup
  normalize = True
  response = "outcome"

  headers = data.pop(0)
  resp_index = headers.index(response)

  if normalize:
    cols = [[row[i] for row in data] for i, h in enumerate(headers)]
    cols = [norm(col) if (all(map(is_num, col)) and i!=resp_index) else col
                    for i, col in enumerate(cols)]
    data = [[col[i] for col in cols] for i in xrange(len(data))]

  return [headers]+data


def _default_postprocessor(splits, headers):
  """
  The default post-processor to be used in model-generation
  after the data-splitting process. Removes unnecessary headers
  and sets values to those expected by the model.
  """

  # Remove unnecessary fields from the data.
  blacklist_fields = [
  "XXXtitle", "XXXinorg1", "XXXinorg2", "XXXinorg3",
  "XXXorg1", "XXXorg2", "XXXoxlike1"
  ]
  blacklist = {headers.index(field)
                 for field in blacklist_fields if field in headers}

  headers = [field for i, field in enumerate(headers) if i not in blacklist]
  splits = { key:[ [e for i,e in enumerate(row) if i not in blacklist]
                        for row in data]
                          for key, data in splits.items()}


  # Convert any stringy-bools in the data to integers.
  bool_map = {
  "yes":1,
  "no":0,
  "?":-1
  }


  for key, data in splits.items():
    for i, row in enumerate(data):
      for j, elem in enumerate(row):
        if elem in bool_map:
          splits[key][i][j] = bool_map[elem]


  return (splits, headers)


def gen_model(title, description, data=None, debug=False, active=False):

  from DRP.models import ModelStats
  from DRP.retrievalFunctions import get_valid_data
  from DRP.model_building.rxn_calculator import headers

  # Prepare the default data if it is unavailable.
  if data is None:
    if debug:
      print "Gathering default data..."
    data = get_valid_data()
    data = [headers]+[d.get_calculations_list() for d in data]


  splitter = None # Enables the "default" splitter.
  #from DRP.model_building.splitters import naive_split as splitter

  model = ModelStats()
  model.construct(title, data,
                  description=description,
                  active=active,
                  preprocessor=_default_preprocessor,
                  postprocessor=_default_postprocessor,
                  splitter=splitter,
                  tool="svc",
                  library="weka",
                  debug=debug)
  model.summary()

  return model



def build_previous_model(model_name, model_description, date, data=None, active=False):
  """
  Constructs a model from the data available on a given date.
  """

  from DRP.model_building import model_methods
  from DRP.retrievalFunctions import filter_by_date, filter_existing_calcs
  from DRP.models import Data

  if not data:
    data = filter_existing_calcs(Data.objects.all())

  filtered = filter_by_date(data, date, "previous")
  model_methods.gen_model(model_name, model_description, data=filtered, active=active)


def retrogenerateModel(date):
  """
  A convenient wrapper to generate a model using data available on a given date.
  """

  title = "Retrogenerated_{}".format(date.replace(" ","_").replace("-","_"))
  description = "A model retrogenerated from the data available on {}".format(date)
  build_previous_model(title, description, date)


def retrogenerateModels():
  """
  Constructs a Learning Curve based on time by repeatedly retrogenerating models.
  """

  def dateRange(start, interval):
    import datetime
    from dateutil.relativedelta import relativedelta

    if interval == "months":
      interval = relativedelta(months=1)
    elif interval == "days":
      interval = relativedelta(days=10)
    else:
      raise Exception("No interval prepared for '{}'".format(interval))

    end = datetime.datetime.now()
    current = start
    while current < end:
      yield current
      current += interval
    yield end # Finally, use *all* the data up to the current time.

  from DRP.models import Data
  earliest_datum = Data.objects.order_by("creation_time_dt")[0]
  start = earliest_datum.creation_time_dt

  for date in dateRange(start, "months"):
    date_string = date.strftime("%m-%d-%Y")
    print "Retrogenerating model from {}".format(date_string)
    retrogenerateModel(date_string)


if __name__=="__main__":
  retrogenerateModels()
