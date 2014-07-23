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

path_to_vis_data_file = "/home/ntien/workspace/dark-reaction-site/DRP/views/vis_data.json" 

@login_required
def get_graph_data(request):
  import os.path 

  #If vis_data is already created and up to date, just return that file 
  if os.path.exists(path_to_vis_data_file) == True:
    print "vis_data exists"
    file_data = open(path_to_vis_data_file)
    deserialized_data = json.load(file_data)
    json_formatted_to_string = json.dumps(deserialized_data) 
    file_data.close() 
    return HttpResponse(json_formatted_to_string, content_type="application/json")   
  
  #If vis_data is not created or not up to date, write new vis_data with current data and return that
  else:
    print "vis_data does not exist" 
    #Only grab reactions that have DataCalc objects already generated
    data = Data.objects.filter(~Q(calculations=None))[:1000] 
	    
    #Grab all data id  
    did = [datum.id for datum in data]  
    print "just  finished querying for data objects and appending dids" 	  
    #Append the data id of each Reaction(DataCalc object) to the end of the row (will be the last field)
    #Works because data ids' and expanded_data's reactions are in the same order.  
    expanded_data = expand_data(data)  
    for i in xrange(len(expanded_data)):
      expanded_data[i] = expanded_data[i] + [did[i]]   
    #Send headers and expanded data in rough CSV form to vis/dataMatrix function (cleans data) in preparation to be put into JSON object form 
    headers = get_expanded_headers() + ["id"]   
    cleaned_matrix = dataMatrix([headers] + expanded_data)
    #Send cleaned matrix to dataMatrixToGraph.myGraph (puts into correct "object" format for  d3 graph)	  
    print "Just finished creating dataMatrix to graph" 
    matrix_formatted_for_vis = myGraph(cleaned_matrix)
    matrix_prepped_for_json = matrix_formatted_for_vis.writeJson()
    #this_file_path = os.path.dirname(os.path.realpath((__file__))) 
    #completeName = os.path.join(this_file_path, "vis_data")   
    print "Just finished creating myGraph object" 
	  
    vis_file = create_vis_data_file(matrix_prepped_for_json) 
    print "just finished creating vis_data file"  
    file_data = open(path_to_vis_data_file) 
    deserialized_data = json.load(file_data)
    json_formatted_to_string = json.dumps(deserialized_data) 
    file_data.close() 
    return HttpResponse(json_formatted_to_string, content_type="application/json")   
  
  #return HttpResponse(json.dumps(matrix_to_json), content_type="application/json")

def create_vis_data_file(data_to_file):
  print path_to_vis_data_file
  with open(path_to_vis_data_file, "w") as outfile:
    dump = json.dump(data_to_file, outfile) 
  return dump 

