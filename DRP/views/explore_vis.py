from django.http import HttpResponse
from django.shortcuts import render
from django.views.static import serve

from DRP.settings import BASE_DIR

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

path_to_vis_data_file = BASE_DIR + "/DRP/views/vis_data.json" 

@login_required
def get_graph_data(request):
  import os.path 

  path_to_nodePositions = BASE_DIR + "/DRP/vis/nodePositions.json"
  path_to_links = BASE_DIR + "/DRP/vis/linkIndices.json" 
  path_to_label_dict = BASE_DIR + "/DRP/vis/label_node_dict.json" 
  #If vis_data is already created and up to date, just return that file 
  if os.path.exists(path_to_nodePositions) and os.path.exists(path_to_links): 
    with open(path_to_nodePositions, "r") as f:
      nodes = json.load(f) 
   
    with open(path_to_links, "r") as f:
      links = json.load(f)

    if os.path.exists(path_to_label_dict):
      with open(path_to_label_dict, "r") as f:
        node_clusters = json.load(f) 
    else: 
      node_clusters = create_node_clusters_for_labels(links, nodes)
      with open(path_to_label_dict, "w") as outfile:
        json.dump(node_clusters, outfile)       
    print node_clusters[0]  
    response = {"nodes": nodes, "links": links, "skipTicks": "True"} 
    return HttpResponse(json.dumps(response), content_type="application/json")   
  
  #If vis_data is not created or not up to date, write new vis_data with current data and return that
  elif os.path.exists(path_to_vis_data_file):
    
    print "vis_data exists" 
    with open(path_to_vis_data_file, "r") as f:
      deserialized_data = json.load(f)
       
    node_clusters = create_node_clusters_for_labels(deserialized_data["links"], deserialized_data["nodes"]) 

    return HttpResponse(json.dumps(deserialized_data), content_type="application/json") 
  else:
    print "vis_data does not exist" 
    #Only grab reactions that have DataCalc objects already generated
    data = Data.objects.filter(~Q(calculations=None)) 
	    
    #Grab all data ids  
    dids= [datum.id for datum in data]  
    print "just  finished querying for data objects and appending dids" 	  
    #Append the data id of each Reaction(DataCalc object) to the end of the row (will be the last field)
    #Works because data ids' and expanded_data's reactions are in the same order.  
    expanded_data = expand_data(data)  
    for i in xrange(len(expanded_data)):
      expanded_data[i] = expanded_data[i] + [dids[i]]   
    #Send headers and expanded data in rough CSV form to datamatrix.py.dataMatrix function (cleans data) in preparation to be put into JSON object form 
    headers = get_expanded_headers() + ["id"]   
    cleaned_matrix = dataMatrix([headers] + expanded_data)
    cleaned_matrix.removeCorrelatedLinregs(0,0)  
    #Send cleaned matrix to dataMatrixToGraph.myGraph (puts into correct "object" format for  d3 graph)	  
    print "Just finished creating dataMatrix to graph" 
    matrix_formatted_for_vis = myGraph(cleaned_matrix)
    matrix_prepped_for_json = matrix_formatted_for_vis.writeJson()
    #this_file_path = os.path.dirname(os.path.realpath((__file__))) 
    #completeName = os.path.join(this_file_path, "vis_data")   
    print "Just finished creating myGraph object" 
	  
    vis_file = create_vis_data_file(matrix_prepped_for_json) 
    print "just finished creating vis_data file"  
        
    with open(path_to_vis_data_file, "r") as f:
        deserialized_data = json.load(f) 

    nodes_clusters= create_node_clusters_for_labels(deserialized_data["links"], deserialized_data["nodes"]) 

    return HttpResponse(json.dumps(deserialized_data), content_type="application/json")   
 

def create_vis_data_file(data_to_file):
  with open(path_to_vis_data_file, "w") as outfile:
    dump = json.dump(data_to_file, outfile) 
  return dump 


def create_node_clusters_for_labels(links, nodes):
  #first create a dictionary of nodes with ids, links (source, target values) and pageranks
  #because we need to sort the links(target, source values) by their associated pagerank and also
  #have access to the link's associated node id (links is in the same order as nodes, as they are both  #a list of dictionaries/objects, but must be combined into one list of objects in order to 
  # math pagerank and id with link information. 
  nodes_dict = [] 
  for i in range(len(nodes)):
    nodes_dict.append({
     "id":nodes[i]["id"],
     "source": links[i]["source"],
     "target": links[i]["target"], 
     "pagerank": nodes[i]["pagerank"],
     "inorg1": nodes[i]["inorg1"],
     "inorg2": nodes[i]["inorg2"] 
     })
  from operator import itemgetter 
  #now the dict has been created, sort it by pagerank (denotes the most "common" reactions/the ones
  #with the most links that should therefore be in the center for labelling purposes
  sorted_dict = sorted(nodes_dict, key=itemgetter("pagerank"), reverse=True)
  #now create the clusters of nodes that are linked directly (
  clusters = find_node_clusters(sorted_dict)
  filtered_clusters = filter_label_clusters(clusters)    
  return filtered_clusters

def find_node_clusters(dictionary):
  neighbors = [] 
  while len(dictionary) > 1:
    for element in dictionary:
      cluster = []
      centerNode = dictionary.pop(0)
      for item in dictionary:
        if check_inorgs(centerNode, item) == True: 
          if item["target"] == centerNode["source"]:
	    cluster.append(item)
	    dictionary.pop(dictionary.index(item)) 
	  else:
	    for entry in cluster:
              if item["target"] == entry["source"]: 
                cluster.append(item)
	        index = dictionary.index(item) 
	        dictionary.pop(index)
	        item = dictionary[index]
      for entry in cluster:
        for item in dictionary:
   	  if item["target"] == entry["source"] and check_inorgs(entry, item) == True:
	    cluster.append(item)
	    dictionary.remove(item)
      cluster.append(centerNode) 
      neighbors.append(cluster)
  else:
   return neighbors  

def filter_label_clusters(cluster):
  filtered_cluster = [item for item in cluster if len(item) > 1] 
  return filtered_cluster
 
#    neighbors.append([item for item in dictionary if item["target"] == centerNode["source"] and check_inorgs(centerNode, item, full_nodes) == True]) 
    #remove these neighbors from the dictionary

#This should find all the neighbors (nodes with direct links (targets match main node's source) and same inorg compounds of the main node, and then all the neighbors of each neighbor (only stopping when there are no more neighbors (and the length of neighbors should continue growing until all neighbors have been found)  

def check_inorgs(mainNode, neighbor):
  if mainNode["inorg1"] == neighbor["inorg1"] and mainNode["inorg2"] == neighbor["inorg2"]:
    return True
  else:
    return False 

def store_graph(request):
  nodeData = json.loads(request.POST["nodes"])
  linkData = json.loads(request.POST["links"]) 
  path_to_nodePositions = BASE_DIR + "/DRP/vis/nodePositions.json"
  path_to_links = BASE_DIR + "/DRP/vis/linkIndices.json" 
  
  with open(path_to_nodePositions, "w") as f:
    json.dump(nodeData, f) 
 
  with open(path_to_links, "w") as f:
    json.dump(linkData, f) 
  
  return HttpResponse("OkeyDokey") 


