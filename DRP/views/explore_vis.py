from django.http import HttpResponse
from django.views.decorators.http import require_http_methods
from django.contrib.auth.decorators import login_required
from django.shortcuts import render

from DRP.models import * 
from DRP.retrievalFunctions import * 
from DRP.jsonViews import * 

import csv, json 

# This view should get the data requested by the user (either the lab group's data, public data, or
# both), use jsonViews to convert it to a CSV/json (in the right format for the d3 visualization), and# then render the template with that data. The template, in turn, will render the javascript (d3 vis)
# with the CSV/json file passed to it (the template). 
@login_required

def select_for_vis(request) #This will select the data particular to the lab group of the user 				 	  #(who is requesting the vis), and should also allow the user to add
			    #public data to the data to be rendered  
  u = request.user
  lab_group = u.get_profile().lab_group

#Get the info from the POST request.
  try:
    filters = request.POST.get("filters")
    if filters:
      filters = json.loads(filters)
    model = request.POST.get("model")
    assert model in {"Recommendation", "Data", "Saved", "CompoundEntry"}
  except Exception as e:
    return HttpResponse("Visualization request failed!")

# Now send this 'model' to the jsonViews to format the CSV correctly 

