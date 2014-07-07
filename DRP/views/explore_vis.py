from django.http import HttpResponse
from django.views.decorators.http import require_http_methods
from django.contrib.auth.decorators import login_required
from django.shortcuts import render
from django.views.static import serve

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
  expanded_data = [row[19:-2] + [max(1, row[-1])]  for row in expand_data(data)]#Not all of the data in the row is necessary, hence the [19:-2]  
  print expanded_data 
  headers = get_expanded_headers()
  new_file = [headers] + expanded_data 
  cleaned_matrix = dataMatrix([headers] + expanded_data) 
  matrix_formatted_for_vis = myGraph(cleaned_matrix) 
  matrix_to_json = writeJson(matrix_formatted_for_vis) 
   
#  return HttpResponse(json.dumps(matrix_to_json), mimetype="application/json")
  response = HttpResponse(content_type='text/csv')
  response['Content-Disposition'] = 'attachment; filename="somefilename.csv" '  
  result = csv.writer(response)
  for row in new_file:
    result.writerow(row)
  return result 
