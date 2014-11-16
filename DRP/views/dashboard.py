## # # # # # # # # # # # # # # # # #
 # # # #  Dashboard Views  # # # # #
# # # # # # # # # # # # # # # # # # #

#Necessary Imports:
from django.http import HttpResponse
from django.shortcuts import render
from django.contrib.auth.decorators import login_required
import json, datetime

from DRP.retrievalFunctions import *
from DRP.models import *


"""
#LEGACY: Used to return dashboard info based on time rather than on version.
#Helper function for get_fields_as_json.
def get_field_tuple(stat, entry):
  #Note: in order to get milliseconds since the epoch, we need a TimeDelta object.
  seconds = int((entry.datetime - datetime.datetime(1970,1,1)).total_seconds()*1000)
  value = getattr(entry, stat)
  return [seconds, value]
"""

def get_fields_as_json(models):
  #Variable Setup.

  # The D3 library needs the following format:
  # [ { "key": "Label for these datums.",
  #     "values": [ [x-value1, y-value1], [x-value2, y-value2], ... ]
  #   },
  #   { ... }
  # ]

  stat_counter = {}
  for i, model in enumerate(models):

    try:
      stats =  model.stats()
    except Exception as e:
      # If the model isn't loadable, mark it as not loadable.
      print "{} was not usable: {}".format(model, e)
      model.check_usability()
      continue

    for stat, val in stats.items():
      if stat not in stat_counter:
        stat_counter[stat] = {"key":stat, "values":[]}

      stat_counter[stat]["values"].append([i, val])

  results = [key_vals for stat, key_vals in stat_counter.items()]
  return results

def get_stats_json(request):
  #Grab all of the model_stats.
  from DRP.retrievalFunctions import get_usable_models
  models = get_usable_models()

  #Convert the data into a JSON format.
  lines = get_fields_as_json(models)
  labels = [model.description for model in models]
  data = json.dumps({"lines":lines, "model_labels":labels})

  #Send the JSON back to the client.
  return HttpResponse(data, mimetype="application/json")

@login_required
def get_dashboard(request):
  return render(request, 'global_page.html', {
   "template": "dashboard",
  })
