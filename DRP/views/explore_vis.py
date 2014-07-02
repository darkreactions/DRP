from django.http import HttpResponse
from django.views.decorators.http import require_http_methods
from django.contrib.auth.decorators import login_required
from django.shortcuts import render

from DRP.models import * 
from DRP.retrievalFunctions import * 
from DRP.views.vis.clustering import *
from DRP.views.vis.datamatrix import *
from DRP.views.vis.datamatrixToGraph import *
from DRP. views.vis.get_data import *
from DRP.views.vis.kdtree import * 

import csv, json 

# This view should get the data requested by the user (either the lab group's data, public data, or
# both), use jsonViews to convert it to a CSV/json (in the right format for the d3 visualization), and# then render the template with that data. The template, in turn, will render the javascript (d3 vis)
# with the CSV/json file passed to it (the template). 
@login_required

def select_for_vis(request) #This will select the data particular to the lab group of the user 				 	  #(who is requesting the vis), and should also allow the user to add
			    #public data to the data to be rendered  
  u = request.user
  lab_group = u.get_profile().lab_group

# result = json.dumps(data)
# context = RequestContext(request, {
#    'result': result,
#})
# return HttpResponse(template.render(context)) 



# Now send this 'model' to the jsonViews to format the CSV correctly 
  csv = write_expanded_data_to_csv() 
  cleaned_matrix = datamatrix.dataMatrix(csv)
  matrix_formatted_for_vis = datamatrixToGraph.myGraph(cleaned_matrix) 
  matrix_to_json = datamatrixToGraph.writeJson(matrix_formatted_for_vis) 
  context = RequestContext(request, {
    'matrix_to_json':matrix_to_json,
  }) 
  return HttpResponse(template.render(context)) 
     
