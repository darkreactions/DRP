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

  """
  # Add bad recommendations to the dataset as 0-outcome reactions.
  from DRP.models import Recommendation
  recs = Recommendation.objects.filter(nonsense=True)
  clean_recs = []
  for i, rec in enumerate(recs):

    try:
      rec.get_calculations_list()
      clean_recs.append(rec)
    except Exception as e:
      pass

  print "NUM CLEAN RECS: {} / {}".format(len(clean_recs), len(recs))
  data += clean_recs
  """


  return data

def generate_avg(title, data, iterations=3, only_keep_avg=True, construct_kwargs={}):

  from DRP.model_building.confusion_table import get_avg_confusion_dict
  from DRP.models import ModelStats

  debug = "debug" in construct_kwargs and construct_kwargs["debug"]==True

  if debug:
    print "Performing {} model-gen iterations...".format(iterations)

  # Construct multiple `iterations` of models
  model_stats = []
  for i in xrange(iterations):
    model_title = "{}_{}".format(title, i)
    new_model = ModelStats()
    new_model.construct(model_title, data, **construct_kwargs)
    model_stats.append(new_model)
    if debug:
      print ""


  best_model = max(model_stats, key=lambda model: model.accuracy())

  if debug:
    print "Preparing average model..."

  avg_model = ModelStats()
  avg_model.title = title
  avg_model.iterations=iterations

  # Copy some stats directly from the `best_model`.
  copy_from_best = ["headers", "correct_vals", "description", "tags",
                    "filename", "active", "usable", "library", "tool",
                    "response"]
  for field in copy_from_best:
    value = getattr(best_model, field)
    setattr(avg_model, field, value)

  avg_model.tags += " averaged"

  # Set the start and end times of this model as the time taken
  # for the entire sequence of iterations to complete.
  avg_model.start_time = model_stats[0].start_time
  avg_model.end_time = model_stats[-1].end_time

  train_cm = get_avg_confusion_dict(model_stats, table="train")
  avg_model.set_confusion_table(train_cm, table="train")

  test_cm = get_avg_confusion_dict(model_stats, table="test")
  avg_model.set_confusion_table(test_cm, table="test")

  avg_model.save()

  if only_keep_avg:
    for model in model_stats:
      model.active = False
      model.usable = False
      model.save()

  return avg_model



def gen_model(title, description, data=None, test_set=None, force=False,
                                  debug=False, pipeline_test=False,
                                  active=False, tags=""):
  """
  Note: `test_set` is an optional list of data entries that are assumed
        to be pre-processed. If it is set, a test set will not be sampled
        from `data` further on in the pipeline.
  """

  from DRP.retrievalFunctions import get_valid_data
  from DRP.model_building.rxn_calculator import headers
  import random

  # Prepare the default data if it is unavailable.
  if data is None:
    if debug:
      print "Gathering default data..."
    data = list(get_valid_data())

    # Make sure you remark on the filter you're using in the description!
    data = [headers] + research_data_filter(data)

    if debug:
      print "Found {} data...".format(len(data)-1)


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
                  "test_set":test_set,
                  "tool":"svc",
                  "library":"weka",
                  "force":force,
                  "debug":debug,
                  }

  if pipeline_test:
    headers = data.pop(0)
    data = [headers] + random.sample(data, len(data)/50)
    model = generate_avg(title, data, construct_kwargs=construct_kwargs,
                         iterations=1)
  else:
    model = generate_avg(title, data, construct_kwargs=construct_kwargs)

  if debug:
    print "Average model produced:"
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
  from DRP.preprocessors import default_preprocessor
  from DRP.model_building.splitters import default_splitter
  import math, datetime, random

  if debug:
    print "Starting at {}".format(datetime.datetime.now().time())

  if data is None:
    if debug:
        print "Gathering default data..."

    # Prepare the default data if it is unavailable.
    data = list(get_valid_data())

    # Make sure you remark on the filter you're using in the description!
    data = [headers] + research_data_filter(data)

    # Take the test-set out of the data (to set aside for the learning curve).
    data = default_preprocessor(data)
    splits = default_splitter(data, headers=headers)
    data = splits["train"]
    test_set = splits["test"]

  else:
    data = list(data)
    test_set = None

  headers = data.pop(0)

  # Randomize the data so that any pre-existing order is not a variable
  #  in experimentation.
  random.shuffle(data)

  bucket_size = int(len(data)*step)
  num_buckets = int(math.ceil(len(data)/float(bucket_size)))

  all_buckets = [data[x:x+bucket_size] for x in xrange(0,len(data),bucket_size)]

  for i in xrange(1, num_buckets+1):

    model_name = "{} ({} of {})".format(name, i, num_buckets)
    model_tag = "learning_curve {} {}".format(curve_tag, i)

    # Grab a random sampling of the data to use.
    union = [bucket for buckets in all_buckets[:i] for bucket in buckets]
    iteration_data = [headers] + union

    if debug:
      print "Generating: \"{}\" (size: {})".format(model_name, len(union))

    # Generate the model.
    gen_model(model_name, description, tags=model_tag,
                                       force=force,
                                       data=iteration_data,
                                       test_set=test_set,
                                       debug=gen_debug)

    if debug:
      print "\tDone: {}\n".format(datetime.datetime.now().time())



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
