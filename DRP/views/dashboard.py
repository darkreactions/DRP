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

def get_fields_as_json(model_stats):
  #Variable Setup.

  # The D3 library needs the following format:
  # [ { "key": "Label for these datums.",
  #     "values": [ [x-value1, y-value1], [x-value2, y-value2], ... ]
  #   },
  #   { ... }
  # ]

  results = [
              {
                "key":"Test Accuracy",
                "values":[[i, model.test_accuracy()] for i, model in enumerate(model_stats)],
              },
              {
                "key":"Test Precision",
                "values":[[i, model.test_precision()] for i, model in enumerate(model_stats)],
              },
              {
                "key":"False Positive Rate",
                "values":[[i, model.rate("false_positive")] for i, model in enumerate(model_stats)],
              },
              {
                "key":"False Negative Rate",
                "values":[[i, model.rate("false_negative")] for i, model in enumerate(model_stats)],
              },
              {
                "key":"User Satisfaction",
                "values":[[i, model.user_satisfaction()] for i, model in enumerate(model_stats)],
              },
            ]

  return results

def get_model_array(model_stats):
  return [entry.description for entry in model_stats]

def get_stats_json(request):
  #Grab all of the model_stats.
  model_stats = ModelStats.objects.all().order_by("datetime")

  #Convert the data into a JSON format.
  lines = get_fields_as_json(model_stats)
  model_labels = get_model_array(model_stats)
  data = json.dumps({"lines":lines, "model_labels":model_labels})

  #Send the JSON back to the client.
  return HttpResponse(data, mimetype="application/json")

@login_required
def get_dashboard(request):
  return render(request, 'global_page.html', {
   "template": "dashboard",
  })
