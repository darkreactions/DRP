#!/usr/bin/env python

#Grab the Django settings if they aren't already set.
import os, sys
full_path = os.path.dirname(os.path.realpath(__file__))+"/"
django_path = full_path[:full_path.rfind("/DRP/")]
if django_path not in sys.path:
  sys.path = [django_path] + sys.path
  os.environ['DJANGO_SETTINGS_MODULE'] = 'DRP.settings'


def build_previous_model(model_name, model_description, date, data=None):
  """
  Constructs a model from the data available on a given date.
  """

  from DRP.model_building import model_methods
  from DRP.retrievalFunctions import filter_by_date, filter_existing_calcs
  from DRP.models import Data

  if not data:
    data = filter_existing_calcs(Data.objects.all())

  filtered = filter_by_date(data, date, "previous")
  model_methods.gen_model(model_name, model_description, data=filtered)

def retrogenerateModel(date):
  """
  A convenient wrapper to generate a model using data available on a given date.
  """

  title = "Retrogenerated {}".format(date)
  description = "A model retrogenerated from the data available on {}".format(date)
  build_previous_model(title, description, date)



def retrogenerateModels():
  """
  Constructs a Learning Curve based on time.
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
 

