# # # # # # # # # # # # # # # # # # #
# # # #  Dashboard Views  # # # # #
# # # # # # # # # # # # # # # # # # #

# Necessary Imports:
from django.http import HttpResponse
from django.shortcuts import render
from django.contrib.auth.decorators import login_required


def get_fields_as_json(models, category="2-test"):
    # Variable Setup.

    # The D3 library needs the following format:
    # [ { "key": "Label for these datums.",
    #     "values": [ [x-value1, y-value1], [x-value2, y-value2], ... ]
    #   },
    #   { ... }
    # ]

    stat_counter = {}

    if not models:
        # If no models are specified, return an empty list of values.
        return []

    for i, model in enumerate(models):

        try:
            stats = model.stats(category=category)
        except Exception as e:
            # If the model isn't loadable, mark it as not loadable.
            print "{} was not usable: {}".format(model, e)
            model.check_usability()
            continue

        for stat, val in stats.items():
            if stat not in stat_counter:
                stat_counter[stat] = {"key": stat, "values": []}

            stat_counter[stat]["values"].append([i, val])

    return stat_counter.values()


@login_required
def get_class_stats_json(request, category="2-test"):

    from DRP.retrievalFunctions import get_usable_models
    from DRP.filters import apply_filters
    import json

    models = get_usable_models()
    models = apply_filters(request, models, model="ModelStats")

    # Convert the data into a JSON format.
    raw = {
        "lines": get_fields_as_json(models, category=category),
        "descriptions": [model.description for model in models],
        "titles": [model.title for model in models],
        "ids": [model.id for model in models],
        "confusionTables": [model.load_confusion_table(normalize=False) for model in models]
    }

    data = json.dumps(raw)

    # Send the JSON back to the client.
    return HttpResponse(data, mimetype="application/json")


@login_required
def get_dashboard(request):

    # Pass the query through the template into the JS.
    current_query = "?" + request.GET.urlencode()

    return render(request, 'global_page.html', {
        "template": "dashboard",
        "current_query": current_query,
    })
