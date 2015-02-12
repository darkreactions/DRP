#!/usr/bin/env python

#Grab the Django settings if they aren't already set.
import os, sys
full_path = os.path.dirname(os.path.realpath(__file__))+"/"
django_path = full_path[:full_path.rfind("/DRP/")]
if django_path not in sys.path:
  sys.path = [django_path] + sys.path
  os.environ['DJANGO_SETTINGS_MODULE'] = 'DRP.settings'

def research_data_filter(data):
  """
  Allows easily filtering operations on the data.
  Use for research but not in the standard model pipeline.

  Make sure you take note of the purpose of this filter
  in the model description!
  """


  # Developers: Put any processing steps here.


  return data


def gen_model(title, description, data=None, debug=False, active=False, tags=""):

  from DRP.models import ModelStats
  from DRP.retrievalFunctions import get_valid_data
  from DRP.model_building.rxn_calculator import headers

  # Prepare the default data if it is unavailable.
  if data is None:
    if debug:
      print "Gathering default data..."
    data = get_valid_data()

    # Make sure you remark on the filter you're using in the description!
    data = research_data_filter(data)

    data = [headers]+[d.get_calculations_list() for d in data]


  # If `splitter` is set to `None`, the default splitter will be used.
  from DRP.preprocessors import default_preprocessor as preprocessor
  from DRP.postprocessors import default_postprocessor as postprocessor
  from DRP.model_building.splitters import default_splitter as splitter

  model = ModelStats()
  model.construct(title, data,
                  description=description,
                  tags=tags,
                  active=active,
                  preprocessor=preprocessor,
                  postprocessor=postprocessor,
                  splitter=splitter,
                  tool="svc",
                  library="weka",
                  debug=debug)
  model.summary()

  return model



def build_model_from_date(model_name, model_description, date, batch_tag, data=None, active=False):
  """
  Constructs a model from the data available on a given date.
  """

  from DRP.retrievalFunctions import filter_by_date, filter_existing_calcs
  from DRP.models import Data

  if not data:
    data = filter_existing_calcs(Data.objects.all())

  filtered = filter_by_date(data, date, "previous")
  tags = "retrogenerated {}".format(batch_tag)

  gen_model(model_name, model_description, data=filtered,
                        active=active, tags=tags)




def retrogenerateModels():
  """
  Constructs a Learning Curve over time by repeatedly retrogenerating models.
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
  import time

  # Create a tag for this batch of recommendations.
  batch_tag = str(int(time.time.now()))

  # Get the epoch datum.
  earliest_datum = Data.objects.order_by("creation_time_dt")[0]
  start = earliest_datum.creation_time_dt

  for date in dateRange(start, "months"):
    date_string = date.strftime("%m-%d-%Y")

    print "Retrogenerating model from {}".format(date_string)

    title = "Retrogenerated_{}".format(date_string.replace(" ","_").replace("-","_"))
    description = "A model generated using data available on {}".format(date_string)
    build_model_from_date(title, description, date_string, batch_tag)


if __name__=="__main__":
  retrogenerateModels()
