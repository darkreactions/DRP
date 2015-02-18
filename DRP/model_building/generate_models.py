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

def generate_avg(self, title, data, iterations=3, only_keep_avg=True,
                                     construct_kwargs={}):
  def _get_avg_confusion_dict(model_stats):
    """
    Returns an average confusion dict from a list of model_stats
    """

    conf_dicts = [m.load_confusion_dict for m in model_stats]
    avg_dict = {}
    for conf_dict in conf_dicts:
      for guess, actuals in conf_dict.items():

        if not guess in avg_dict:
          avg_dict[guess] = {}

        for actual, occurrences in actuals.items():
          if not actual in avg_dict[guess]:
            avg_dict[guess][actual] = 0

          avg_dict[guess][actual] += occurrences

    num_models = float(len(model_stats))
    avg_dict ={guess:{actual:count/num_models for actual,count in actuals.items()}
                     for guess, actuals in avg_dict.items()}

    return avg_dict

  from DRP.models import ModelStats

  # Construct multiple `iterations` of models
  model_stats = [ModelStats().construct(title+str(i), data, **construct_kwargs)
                   for i in xrange(iterations)]

  best_model = max(model_stats, key=lambda model: model.test_accuracy)

  avg_model = ModelStats()
  avg_model.title = title

  # Copy some stats directly from the `best_model`.
  copy_from_best = ["headers", "correct_vals", "description", "tags",
                    "filename", "active", "usable", "library", "tool",
                    "response"]
  for field in copy_from_best:
    value = getattr(best_model, field)
    setattr(avg_model, field, value)

  # Set the start and end times of this model as the time taken
  # for the entire sequence of iterations to complete.
  avg_model.start_time = model_stats[0].start_time
  avg_model.end_time = model_stats[-1].end_time

  avg_model.set_confusion_table(_get_avg_confusion_dict(model_stats))
  avg_model.save()

  if only_keep_avg:
    for model in model_stats:
      model.delete()

  return avg_model



def gen_model(title, description, data=None, force=False, debug=False,
                                             iterations=1,
                                             active=False, tags=""):

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

  construct_kwargs = {
                  "description":description,
                  "tags":tags,
                  "active":active,
                  "preprocessor":preprocessor,
                  "postprocessor":postprocessor,
                  "splitter":splitter,
                  "tool":"random forest",
                  "library":"sklearn",
                  "force":force,
                  "debug":debug
                  }

  model = generate_avg(title, data, iterations=iterations,
                                    construct_kwargs=construct_kwargs)

  if debug:
    model.summary()

  return model



def learning_curve(name, description, curve_tag, data=None,
                                                 force=True,
                                                 step=0.05,
                                                 gen_debug=False,
                                                 debug=False):
  def curve_generator(total_size, step):
    current = step*total_size
    i = 1
    while total_size > current:

      yield i, int(current)

      current += (step*total_size)
      i += 1

    yield i, int(total_size)

  from DRP.retrievalFunctions import get_valid_data
  from DRP.model_building.rxn_calculator import headers
  import random, math, datetime

  if data is None:
    data = get_valid_data()

    # Prepare the default data if it is unavailable.
    if debug:
        print "Gathering default data..."

    # Make sure you remark on the filter you're using in the description!
    data = research_data_filter(data)

    data = [headers]+[d.get_calculations_list() for d in data]

  else:
    data = list(data)

  headers = data.pop(0)

  total_iters = int(math.ceil(1.0/step))

  for i, sample_size in curve_generator( len(data), step):

    model_name = "{}__{}_of_{}".format(name, i, total_iters)
    model_tag = "learning_curve {} {}".format(curve_tag, i)

    # Grab a random sampling of the data to use.
    iteration_data = [headers] + random.sample(data, sample_size)

    if debug:
      print "Generating: {}".format(model_name),

    # Generate the model.
    gen_model(model_name, description, tags=model_tag,
                                       force=force,
                                       data=iteration_data,
                                       debug=gen_debug)

    if debug:
      print " ({})".format(datetime.datetime.now().time())



def build_model_from_date(model_name, model_description, date, batch_tag, data=None):
  """
  Constructs a model from the data available on a given date.
  """

  from DRP.retrievalFunctions import filter_by_date, filter_existing_calcs
  from DRP.models import Data

  if not data:
    data = filter_existing_calcs(Data.objects.all())

  filtered = filter_by_date(data, date, "previous")
  tags = "retrogenerated {}".format(batch_tag)

  gen_model(model_name, model_description, data=filtered, tags=tags)



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
