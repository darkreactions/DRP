from django.http import HttpResponse
from django.shortcuts import render
from django.views.static import serve

from DRP.models import *
from DRP.retrievalFunctions import *
from DRP.vis.datamatrix import *
from DRP.vis.datamatrixToGraph import *
from DRP.vis.get_data import *
from DRP.vis.kdtree import *
from django.contrib.auth.decorators import login_required

from DRP.database_construction import *
from DRP.forms import * 
from DRP.validation import *

from DRP.data_config import CONFIG

import csv, json, string 

# This view should get the data requested by the user (either the lab group's data, public data, or
# both), use jsonViews to convert it to a CSV/json (in the right format for the d3 visualization), and# then render the template with that data. The template, in turn, will render the javascript (d3 vis)
# with the CSV/json file passed to it (the template).
@login_required
def get_graph_data(request):
  #Only grab reactions that have DataCalc objects already generated
  data = Data.objects.filter(~Q(calculations=None))[:100] 
    
  #Grab all data id  
  did = [datum.id for datum in data]  
  
  #Append the data id of each Reaction(DataCalc object) to the end of the row (will be the last field)
  #Works because data ids' and expanded_data's reactions are in the same order. 
  expanded_data = expand_data(data)  
  for i in xrange(len(expanded_data)):
    expanded_data[i] = expanded_data[i] + [did[i]]   
#Seind headers and expanded data in rough CSV form to vis/dataMatrix function (cleans data) in preparation to be put into JSON object form 
  headers = get_expanded_headers() + ["id"]   
  cleaned_matrix = dataMatrix([headers] + expanded_data)
  #Send cleaned matrix to dataMatrixToGraph.myGraph (puts into correct "object" format for  d3 graph)
  matrix_formatted_for_vis = myGraph(cleaned_matrix)
  matrix_to_json = matrix_formatted_for_vis.writeJson()
  print "Printing path"
  import os.path
  print os.path.isfile("~/DRP/vis/vis_data")    
  this_file_path = os.path.dirname(os.path.realpath((__file__))) 
  completeName = os.path.join(this_file_path, "vis_data")   
  
  #If vis_data is already created and up to date, just return that file 
  if os.path.isfile("/DRP/vis/vis_data") == True:
    print "vis_data exists"
    print completeName 
    return HttpResponse(json.dumps(matrix_to_json), content_type="application/json")   

  #If vis_data is not created or not up to date, write new vis_data with current data and return that
  else:
    if os.path.isfile("vis_data") == False: 
      print "here" 
      print os.path:  
      with open(completeName, "w") as outfile:
        dump = json.dumps(matrix_to_json, outfile, indent = 2) #, ensure_ascii=False)   
      return HttpResponse(dump, content_type="application/json")   
  
  #return HttpResponse(json.dumps(matrix_to_json), content_type="application/json")



