## # # # # # # # # # # # # # # # # #
 # # # #  Dashboard Views  # # # # #
# # # # # # # # # # # # # # # # # # #

#Necessary Imports:
from django.http import HttpResponse
from django.shortcuts import render
from django.contrib.auth.decorators import login_required


def get_fields_as_json(models, classes=4):
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
      stats =  model.stats(classes=classes)
    except Exception as e:
      # If the model isn't loadable, mark it as not loadable.
      print "{} was not usable: {}".format(model, e)
      model.check_usability()
      continue


    for stat, val in stats.items():
      if stat not in stat_counter:
        stat_counter[stat] = {"key":stat, "values":[]}

      stat_counter[stat]["values"].append([i, val])

  size_tups = stat_counter["Test Size"]["values"]
  max_size = float(max([size for index, size in size_tups]))
  new_size_tups = [(index, size/max_size) for index, size in size_tups]
  stat_counter["Test Size"]["values"] = new_size_tups


  return stat_counter.values()


def get_class_stats_json(request, classes=4):

  from DRP.retrievalFunctions import get_usable_models
  import json

  models = get_usable_models()

  #Convert the data into a JSON format.
  raw = {
    "lines":get_fields_as_json(models, classes=classes),
    "descriptions":[model.description for model in models],
    "titles":[model.title for model in models],
    "confusionTables":[model.load_confusion_table(normalize=False) for model in models]
  }

  data = json.dumps(raw)

  #Send the JSON back to the client.
  return HttpResponse(data, mimetype="application/json")


@login_required
def get_dashboard(request):
  return render(request, 'global_page.html', {
   "template": "dashboard",
  })
