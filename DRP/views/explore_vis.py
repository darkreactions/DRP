from django.http import HttpResponse
from django.views.decorators.http import require_http_methods
from django.contrib.auth.decorators import login_required
from django.shortcuts import render

from DRP.models import * 
from DRP.retrievalFunctions import * 
from DRP.vis.datamatrix import *
from DRP.vis.datamatrixToGraph import *
from DRP.vis.get_data import *
from DRP.vis.kdtree import * 

import csv, json 

# This view should get the data requested by the user (either the lab group's data, public data, or
# both), use jsonViews to convert it to a CSV/json (in the right format for the d3 visualization), and# then render the template with that data. The template, in turn, will render the javascript (d3 vis)
# with the CSV/json file passed to it (the template). 
@login_required
def get_graph_data(request): 
  data = Data.objects.filter(~Q(calculations=None)) #Only grab reactions that have DataCalc objects already generated
  expanded_data = [row[19:-2] + [max(1, row[-1])]  for row in expand_data(data)] 
  print expanded_data 
  headers = get_expanded_headers() 
  cleaned_matrix = dataMatrix([headers] + expanded_data) 
  matrix_formatted_for_vis = myGraph(cleaned_matrix) 
  matrix_to_json = writeJson(matrix_formatted_for_vis) 
   
  return HttpResponse(json.dumps(matrix_to_json), mimetype="application/json")




def get_vis(request): #This will select the data particular to the lab group of the user 				 	  #(who is requesting the vis), and should also allow the user to add
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
     
