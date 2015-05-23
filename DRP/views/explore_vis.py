'''
This view should get the data requested by the user (either the lab group's data, public data, or both), use jsonViews to convert it to a CSV/json (in the right format for the d3 visualization), and then render the template with that    data. The template, in turn, will render the javascript (d3 vis) with the CSV/json file passed to it (the template). Currently, however, a lot of this is hardcoded.

vis_data.json ----> has all of the nodes and links created by the KDtree nearest neighbor file from the data (so new data is not included/graph is not dynamic) because it takes forever to create the nodes/links from the database. Should look at KDtree code to speed it up. The links are the target and sources (the edges between nodes/which nodes are linked). The nodes have all the desired attributes EXCEPT x and y positions (which were stored in a file generated from a run of the D3.js graph in order to avoid needing to recalculate positions every time).

nodePositions.json ---->

linkIndices.json ----->

label_node_dict.json ----->
'''

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

from DRP.compoundGuideFunctions import translate_reactants

import csv, json, string

#Global Variables
FIRST_CLUSTER_ABBREVS = {"Oxovanadium(2+) sulfate": "VOSO4", "sodium metavanadate": "NaVO3", "ammonium metavanadate": "NH4VO3", "vanadium(V) oxide": "V2O5", "molybdenum trioxide": "MoO3", "gallium trinitrate": "Ga(NO3)3", "sodium molybdate": "Na2MoO4", "potassium vanadium trioxide": "Potassium Vanadium Trioxide", "sodium vanadium trioxide": "Sodium Vanadium Trioxide", "Ga2O3": "Ga2O3", "potassium metavanadate": "KO3"}

SECOND_CLUSTER_ABBREVS = {
"Oxovanadium(2+) sulfate Selenium dioxide": "VOSO4 + SeO2",
"sodium tellurite sodium metavanadate": "Na2TeO3 + NaVO3",
"molybdenum trioxide": "MoO3",
"sodium metavanadate": "NaVO3",
"sodium metavanadate sodium tellurite": "NaVO3 + Na2TeO3",
"ammonium metavanadate Selenium dioxide": "NH4VO3 + SeO2",
"sodium metavanadate Selenium dioxide": "NaVO3 + SeO2",
"sodium tellurite ammonium metavanadate": "Na2TeO3 + NH4VO3",
"sodium metavanadate tellurium dioxide": "NaVO3 + TeO2",
"Ga2O3": "Ga2O3",
"vanadium(V) oxide sodium tellurite": "V2O5 + Na2TeO3",
"vanadium(V) oxide Selenium dioxide": "V2O5 + SeO2",
"ZnO": "ZnO","Zn(NO3)2": "Zn(NO3)2",
"gallium trinitrate": "Ga(NO3)3",
"potassium vanadium trioxide Selenium dioxide": "Potassium Vanadium Trioxide + SeO2",
"sodium molybdate": "Na2MoO4","sodium tellurite vanadium(V) oxide": "Na2TeO3 + V2O5",
"sodium vanadium trioxide Selenium dioxide": "Sodium Vanadium Trioxide + SeO2",
"Oxovanadium(2+) sulfate Potassium Dichromate": "VOSO4 + K2Cr2O7",
"sodium fluoride Selenium dioxide": "Sodium Fluoride + SeO2",
"Oxovanadium(2+) sulfate sodium tellurite": "VOSO4 + Na2TeO3",
"vanadium(V) oxide": "V2O5",
"sodium metavanadate Potassium Dichromate": "NaVO3 + K2Cr2O7",
"Oxovanadium(2+) sulfate": "VOSO4",
"ammonium metavanadate": "NH4VO3",
"Selenium dioxide": "SeO3"
}


VIS_DATA_PATH = BASE_DIR + "/DRP/vis/vis_data.json"
NODEPOSITIONS_PATH  = BASE_DIR + "/DRP/vis/nodePositions.json"
NODENOOUTCOMES = BASE_DIR + "/DRP/vis/nodePosNoOutcome.json"
LINKS_PATH = BASE_DIR + "/DRP/vis/new_linkIndices.json"
LABEL_DICT_PATH = BASE_DIR + "/DRP/vis/label_node_dict.json"

COLORS =["#a6cee3", "#1f78b4"," #b2df8a"," #33a02c", "#fb9a99", "#e31a1c", "#fdbf6f", "#b15928", "#cab2d6","#6a3d9a", "#ffff99", "#ff7f00"]

@login_required
def get_graph_data(request):
  import os.path

  #If vis_data is already created and up to date, just return that file
  if os.path.exists(NODENOOUTCOMES) and os.path.exists(LINKS_PATH):
    with open(NODENOOUTCOMES, "r") as f:
      nodes = json.load(f)
    with open(LINKS_PATH, "r") as f:
      links = json.load(f)
    #If the node clusters have already been created, load that file. Otherwise, create node clusters for second and first tier clusters
    if os.path.exists(LABEL_DICT_PATH):
      with open(LABEL_DICT_PATH, "r") as f:
        node_clusters = json.load(f)
    else:
      node_clusters = create_node_clusters_for_labels(links, nodes)
      with open(LABEL_DICT_PATH, "w") as outfile:
        json.dump(node_clusters, outfile)

    #This should definitely exist if nodePositions and links exist/this is just here because nodeElements were acting funky when using nodePositions data in the javascript
  if os.path.exists(VIS_DATA_PATH):
    with open(VIS_DATA_PATH, "r") as f:
      deserialized_data = json.load(f)
    rawNodes = deserialized_data["nodes"]
    #links = deserialized_data["links"]
    if os.path.exists(LINKS_PATH):
      with open(LINKS_PATH, "r") as f:
        links = json.load(f)
    #Here is the code that deals with the randomly missing outcome attributes in nodePositions
    no_outcome = []
    for i in xrange(len(nodes)):
      if "outcome" not in nodes[i]:
        no_outcome.append(nodes[i])

    no_outcome_ids = []
    for i in xrange(len(no_outcome)):
      no_outcome_ids.append(no_outcome[i]["id"])

    raw_outcomes = []
    for i in xrange(len(rawNodes)):
      if rawNodes[i]["id"] in no_outcome_ids:
        raw_outcomes.append(rawNodes[i])

    #Here I am grabbing list of all abbreviated compounds 
     
    from DRP.models import collect_CG_name_pairs
    u = request.user
    lab_group = u.get_profile().lab_group
  
    name_pairs = collect_CG_name_pairs(lab_group)
    name_pairs = { value:key for key,value in name_pairs.items() }    


    with_outcomes = []
    for i in xrange(len(no_outcome)):
      id_ = no_outcome[i]["id"]
      for j in xrange(len(raw_outcomes)):
        if raw_outcomes[j]["id"] == id_:
          outcome_node = raw_outcomes[j]
      with_outcomes.append({
        "purity": outcome_node["purity"],
        "target": no_outcome[i]["target"],
        "color": no_outcome[i]["color"],
        "label1": no_outcome[i]["label1"],
        "label2": no_outcome[i]["label2"],
        "source": no_outcome[i]["source"],
        "color2": no_outcome[i]["color2"],
        "inorg1": no_outcome[i]["inorg1"],
        "inorg1_abbrev": name_pairs[no_outcome[i]["inorg1"]],
        "y": no_outcome[i]["y"],
        "x": no_outcome[i]["x"],
        "pagerank": no_outcome[i]["pagerank"],
        "id": no_outcome[i]["id"],
        "outcome": outcome_node["outcome"],
        "ref": outcome_node["ref"]
        })
      if "inorg2" in no_outcome[i]:
        if no_outcome[i]["inorg2"] != "none":
          with_outcomes.append({
          "inorg2": no_outcome[i]["inorg2"],
          "inorg2_abbrev": name_pairs[no_outcome[i]["inorg2"]]})   
    
    if os.path.exists(LABEL_DICT_PATH):
      with open(LABEL_DICT_PATH, "r") as f:
        clusters_with_source_target = json.load(f)

    with_outcomes[:] = [ x for x in with_outcomes if "id" in x ]
    with_outcomes[:] = [ x for x in with_outcomes if "target" in x ]
    with_outcomes[:] = [ x for x in with_outcomes if "source" in x ]



    nodes = with_outcomes

     #Clusters are the elements that contain all the nodes with matching inorganics (a list of lists of dictionaries)i
    #nodes are the original datum points corresponding to a single reaction( a list of dictionaries)
    #This is for the second tier clusters (clustered by two inorganics)
    clusters = give_positions_to_clusters(nodes, node_clusters)

    #clusters that are objects in a list within a larger list
    clusters_with_colors = assign_colors_to_clusters(COLORS, clusters, "color2")
    final_clusters = make_clusters_into_single_list(clusters_with_colors)

    #Here we are going to make NEW nodes based on final_clusters!!! Necessary because final_cluster are missing nodes that don't fall into any of the clusters, but still needed in vis (must grab x,y coords, inorgs, etc)
    almost_final_nodes = combine_final_clusters_and_nodes(nodes, final_clusters, links) #This seems to be working (color2 appended correctly)
    doublelist = []
    for i in xrange(len(almost_final_nodes)):
      label = almost_final_nodes[i]["label2"]
      if label != "none" and label not in doublelist:
        doublelist.append(label)

    #This is for the first tier clusters (clustered by single inorganic)
    votes = vote_on_inorgs(final_clusters)
    top_inorgs = grab_inorgs(votes)
    firstCluster = first_cluster(almost_final_nodes, top_inorgs)
    first_clusters_with_colors = assign_colors_to_clusters(COLORS, firstCluster, "color")
    final_first_clusters = make_clusters_into_single_list(first_clusters_with_colors)
    final_nodes = give_colors_to_nodes(almost_final_nodes, final_first_clusters) #This also seems to be working as expected, appending color and label1
    singlelist = []
    for i in xrange(len(final_nodes)):
      label = final_nodes[i]["label1"]
      if label != "none" and label not in singlelist:
        singlelist.append(label)

    final_nodes = with_outcomes
    #Here I am replacing the labels with the abbreviated compounds
    for i in xrange(len(final_nodes)):
      label1 = final_nodes[i]["label1"]
      label2 = final_nodes[i]["label2"]
      if label1 in FIRST_CLUSTER_ABBREVS:
        final_nodes[i]["label1"] = FIRST_CLUSTER_ABBREVS[label1]
      if label2 in SECOND_CLUSTER_ABBREVS:
        final_nodes[i]["label2"] = SECOND_CLUSTER_ABBREVS[label2]

    with open(BASE_DIR + "/DRP/vis/completeNodes.json", "w") as outfile:
      dump = json.dump(final_nodes, outfile)
    
    
    clusters = createClusters(final_nodes) 
    clusters1 = clusters[1]
    clusters2 = clusters[0] 
    print len(clusters1)
    labelled_clusters1 = assign_labels_to_all_clusters(node_clusters, clusters1, "label2")
    labelled_clusters2 = assign_labels_to_all_clusters(firstCluster, clusters2, "label1") 

    for i in xrange(len(labelled_clusters1)):
      label = labelled_clusters1[i]["label"]
      if label in SECOND_CLUSTER_ABBREVS:
        labelled_clusters1[i]["label"] = SECOND_CLUSTER_ABBREVS[label]

    for i in xrange(len(labelled_clusters2)):
      label = labelled_clusters2[i]["label"]
      if label in FIRST_CLUSTER_ABBREVS:
        labelled_clusters2[i]["label"] = FIRST_CLUSTER_ABBREVS[label]

    response = {"nodes": final_nodes, "links": links, "skipTicks": "True", "clusters1": labelled_clusters2, "clusters2": labelled_clusters1}

    return HttpResponse(json.dumps(response), content_type="application/json")
   

  else:
    print "vis_data does not exist"
    #Only grab reactions that have DataCalc objects already generated
    all_data = Data.objects.filter(~Q(calculations=None))
    data = all_data[:200] 
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

  return HttpResponse(json.dumps(response), content_type="application/json")


def create_vis_data_file(data_to_file):
  with open(VIS_DATA_PATH, "w") as outfile:
    dump = json.dump(data_to_file, outfile)
  return dump


def create_node_clusters_for_labels(links, nodes):
  #first create a dictionary of nodes with ids, links (source, target values) and pageranks
  #because we need to sort the links(target, source values) by their associated pagerank and also
  #have access to the link's associated node id (links is in the same order as nodes, as they are both a list of dictionaries/objects, but must be combined into one list of objects in order to
  # match pagerank and id with link information).
  nodes_dict = []
  for i in range(len(nodes)):
    nodes_dict.append({
     "id":nodes[i]["id"],
     "source": links[i]["source"],
     "target": links[i]["target"],
     "pagerank": nodes[i]["pagerank"],
     "x": 0,
     "y": 0,
     "color": "none",
     "color2": "none",
     "label1": "none",
     "label2": "none",
     "outcome": nodes[i]["outcome"],
     "ref": nodes[i]["ref"]
     })
    if "inorg1" in nodes[i]:
      nodes_dict[i]["inorg1"] = nodes[i]["inorg1"]
    if "inorg2" in nodes[i]:
      if nodes[i]["inorg2"] != "none":
        nodes_dict[i]["inorg2"] = nodes[i]["inorg2"]
  from operator import itemgetter
  #now the dict has been created, sort it by pagerank (denotes the most "common" reactions/the ones
  #with the most links that should therefore be in the center for labelling purposes
  sorted_dict = sorted(nodes_dict, key=itemgetter("pagerank"), reverse=True)
  #now create the clusters of nodes that are linked directly (
  clusters = find_node_clusters(sorted_dict)
  filtered_clusters = get_rid_of_single_item_clusters(clusters)
  return filtered_clusters

#grab all top-priority inorganics:
def grab_inorgs(list_of_all_inorgs):
  mylist = []
  top_inorgs = []
  for key in list_of_all_inorgs:
    mylist.append(key)
  for i in xrange(len(mylist)):
    if mylist[i] == -1:
      mylist.pop(i)
  for i in xrange(len(mylist)):
    if "va" in mylist[i] or "Va" in mylist[i] or "ga" in mylist[i] or "Ga" in mylist[i] or "mo" in mylist[i] or "Mo" in mylist[i]:
        top_inorgs.append(mylist[i])
  return top_inorgs


#Find all nodes that have the highest priority single organic compounds in common
#Looking at a ranked list of single organic compounds in order to create this first level of the hierarchy
def first_cluster(dictionary_of_all_nodes, singles):
 #singles is the list of all the single, most common/indicative single inorganic compounds
 neighbors = []
 nodedict = dictionary_of_all_nodes
 for i in xrange(len(singles)):
  cluster = []
  for j in xrange(len(nodedict)):
   if nodedict[j]["inorg1"] == singles[i]:
    cluster.append(nodedict[j])
   elif "inorg2" in nodedict[j]:
    if nodedict[j]["inorg2"] == singles[i]:
     cluster.append(nodedict[j])
  neighbors.append(cluster)

 for i in xrange(len(neighbors)):
    if neighbors[i][0]["inorg1"] in singles:
      neighbors[i][0]["label1"] = neighbors[i][0]["inorg1"]
    elif "inorg2" in neighbors[i][0]:
      if neighbors[i][0]["inorg2"] in singles:
       neighbors[i][0]["label1"] = neighbors[i][0]["inorg2"]

 return neighbors

#Finding all nodes that are connected to each other and have both inorgs in common
#Also, finding most common inorgs (single, not in pairs),
def find_node_clusters(dictionary):
  neighbors = []
  while len(dictionary) > 1:
    for element in dictionary:
      cluster = []
      centerNode = dictionary.pop(0)
      if "inorg2" in centerNode and centerNode["inorg2"] != -1.0 and centerNode["inorg2"] != "none": 
        centerNode["label2"] = centerNode["inorg1"] + ", " + centerNode["inorg2"] 
      else:
        centerNode["label2"] = centerNode["inorg1"]
      for item in dictionary:
        #if check_inorgs(centerNode, item) == True:
        if item["target"] == centerNode["source"]:
          cluster.append(item)
          dictionary.pop(dictionary.index(item))
        else:
          for entry in cluster:
            if item["target"] == entry["source"]:
              cluster.append(item)
              index = dictionary.index(item)
              dictionary.pop(index)
              if index < len(dictionary):
                item = dictionary[index] #doing this so the for loop doesn't skip over an element after one is removed
      for entry in cluster:
        for item in dictionary:
   	      if item["target"] == entry["source"] and check_inorgs(entry, item) == True:
	        cluster.append(item)
	        dictionary.remove(item)
      cluster.insert(0,centerNode)
      neighbors.append(cluster)
  else:
   return neighbors

def get_rid_of_single_item_clusters(cluster):
  filtered_cluster = [item for item in cluster if len(item) > 1]
  return filtered_cluster

def make_clusters_into_single_list(clusters):
  result = []
  for i in xrange(len(clusters)):
    for j in xrange(len(clusters[i])):
      result.append(clusters[i][j])
  return result


#This should find all the neighbors (nodes with direct links (targets match main node's source) and same inorg compounds of the main node, and then all the neighbors of each neighbor (only stopping when there are no more neighbors (and the length of neighbors should continue growing until all neighbors have been found)

def check_inorgs(mainNode, neighbor):
  if mainNode["inorg1"] == neighbor["inorg1"]:
    if "inorg2" in mainNode and "inorg2" in neighbor:
      if mainNode["inorg2"] == neighbor["inorg2"]:
        return True
    else:
      return True
  return False

def store_graph(request):
  nodeData = json.loads(request.POST.get("nodes"))
  linkData = json.loads(request.POST.get("links"))
  NODEPOSITIONS_PATH = BASE_DIR + "/DRP/vis/nodePositions.json"
  LINKS_PATH = BASE_DIR + "/DRP/vis/linkIndices.json"

  with open(NODEPOSITIONS_PATH, "w") as f:
    json.dump(nodeData, f)

  with open(LINKS_PATH, "w") as f:
    json.dump(linkData, f)

  return HttpResponse("OkeyDokey")

def give_positions_to_clusters(nodes, clusters):
  for i in xrange(len(clusters)):
    for j in xrange(len(clusters[i])):
      id_num = clusters[i][j]["id"]
      nodeWithPosition = next((x for x in nodes if x["id"] == id_num), None)
      clusters[i][j]["x"] = nodeWithPosition["x"]
      clusters[i][j]["y"] = nodeWithPosition["y"]
  return clusters

#This function is looking at all the clusters for the first level hiearchy (so not all the nodes), and appending a color to the nodes for the nodes that are in the cluster based on 1 single inorganic
def give_colors_to_nodes(nodes, clusters):
    for i in xrange(len(nodes)):
      id_num = nodes[i]["id"]
      colorAndLabel = next((x["color"] + "," + x["label1"] for x in clusters if x["id"] == id_num), None)
      if colorAndLabel != None:
        x = colorAndLabel.index(",")
        nodes[i]["color"] = colorAndLabel[0:x]
        nodes[i]["label1"] = colorAndLabel[x+1:]
    return nodes

#This is doing the same as the above function except for the second level of clusters (based on two inorganics)
def combine_final_clusters_and_nodes(nodes, clusters, links):
    new_nodes = clusters
    for i in range(len(nodes)):
        if nodes[i]["id"] not in [x["id"] for x in clusters]:
            new_node = {
            "id":nodes[i]["id"],
            "source": links[i]["source"], #I am assuming here that links and nodes are in the same order in terms of id numbers **crucial
            "target": links[i]["target"],
            "pagerank": nodes[i]["pagerank"],
            "x": nodes[i]["x"],
            "y": nodes[i]["y"],
            "color": "none",
            "color2": "none",
            "label1": "none",
            "label2": "none",
            "inorg1": nodes[i]["inorg1"], 
            "outcome":nodes[i]["outcome"], 
            "ref": nodes[i]["ref"]
            }
            if "inorg2" in nodes[i] and nodes[i]["inorg2"] != "-1" and nodes[i]["inorg2"] != "none":
              new_node.append({ "inorg2": nodes[i]["inorg2"]}) 
            new_nodes.append(new_node)
    return new_nodes

def reduce_clusters(clusters):
  reduced_clusters = [i[0] for i in clusters]
  return reduced_clusters

def vote_on_inorgs(reduced_clusters):
  inorg_Votes = {}
  for i in xrange(len(reduced_clusters)):
    inorg1 = reduced_clusters[i]["inorg1"]
    if inorg1 in inorg_Votes:
      inorg_Votes[inorg1] += 1
    else:
      inorg_Votes[inorg1] = 1
    if "inorg2" in reduced_clusters[i]:
      inorg2 = reduced_clusters[i]["inorg2"]
      if inorg2 in inorg_Votes:
        inorg_Votes[inorg2] += 1
      else:
        inorg_Votes[inorg2] = 1
  return inorg_Votes

def check_size_of_clusters(clusters):
  small_clusters = 0
  large = 0
  significant = 0
  for i in xrange(len(clusters)):
   if len(clusters[i]) < 15:
    small_clusters += 1
   # print clusters[i][0]["inorg1"]
    if "inorg2" in clusters[i]:
      inorg2 = clusters[i][0]["inorg2"]
   elif len(clusters[i]) > 15 and len(clusters[i]) < 20:
     large += 1
   elif len(clusters[i]) > 20:
     significant += 1
  #print "Large" + str(large)
  #print "Significant" + str(significant)
  return small_clusters


def assign_colors_to_clusters(colors, node_clusters, key):
  for i in xrange(len(node_clusters)):
    for j in xrange(len(node_clusters[i])):
      node_clusters[i][j][key] = colors[i%12]
  return node_clusters

def make_larger_clusters(clusters, total_most_important_inorgs):
  main_clusters = []
  inorgs = total_most_important_inorgs
  centerX = 0
  centerY = 0
  for i in xrange(len(inorgs)):
    centerX = 0
    centerY = 0
    for j in xrange(len(clusters)):
      for k in xrange(len(clusters[j])):
        centerX = (centerX + clusters[j][k]["x"])/2
        centerY = (centerY + clusters[j][k]["y"])/2
        fill = clusters[j][k]["color"]
      if clusters[j][0]["inorg1"] in inorgs[i]:
        main_clusters.append({
          "label": clusters[j][0]["inorg1"],
          "x": centerX,
          "y": centerY,
          "r": 15 * len(clusters[j]),
          "fill": fill
        })
      elif "inorg2" in clusters[j] and clusters[j][0]["inorg2"] in inorgs[i]:
        main_clusters.append({
          "label": clusters[j][0]["inorg2"],
          "x": centerX,
          "y": centerY,
          "r": 15 * len(clusters[j]),
          "fill": fill
          })
  return main_clusters

def make_new_nodePositions():

  with open(VIS_DATA_PATH, "r") as f:
    deserialized_data = json.load(f)
  nodez = deserialized_data["nodes"]
  n = nodes
  print "here"
  new_nodez = []
  for i in xrange(len(n)):
    id_num = n[i]["id"]
    id_node = []
    for j in xrange(len(nodez)):
      if nodez[j]["id"] == id_num:
        id_node = nodez[j]
    new_nodez.append({
    "index": n[i]["index"],
    "target": n[i]["target"],
    "weight": n[i]["weight"],
    "color": n[i]["color"],
    "label1": n[i]["label1"],
    "label2": n[i]["label2"],
    "source": n[i]["source"],
    "color2": n[i]["color2"],
    "y": n[i]["y"],
    "x": n[i]["x"],
    "pagerank": n[i]["pagerank"],
    "px": n[i]["px"],
    "py": n[i]["py"],
    "id": n[i]["id"],
    "outcome": id_node["outcome"],
    "ref": id_node["ref"]
    })
    if "inorg1" in n[i]:
      if n[i]["inorg1"] != "none":
        new_nodez[i]["inorg1"] = n[i]["inorg1"]
    if "inorg2" in n[i]:
      if n[i]["inorg2"] != "none":
        new_nodez[i]["inorg2"] = nodes[i]["inorg2"]
      
  with open(BASE_DIR + "/DRP/vis/nodePosOutcome.json", "w") as outfile:
    json.dump(new_nodez, outfile)
  return HttpResponse("OkeyDokey")

def createClusters(final_nodes):
  nodes = final_nodes
  clusters1 = []
  clusters2 = []
  for i in xrange(len(nodes)):
    if nodes[i]["label1"] != "none": 
      clusters1.append({
        "color": nodes[i]["color"],
        "x": nodes[i]["x"],
        "y": nodes[i]["y"],
        "id": nodes[i]["id"],
        "inorg1": nodes[i]["inorg1"],
        "label": nodes[i]["label1"]  
      })
    if nodes[i]["label2"] != "none": 
      clusters2.append({
        "color": nodes[i]["color2"],
        "x": nodes[i]["x"],
        "y": nodes[i]["y"],
        "id": nodes[i]["id"],
        "inorg1": nodes[i]["inorg1"],
        "label": nodes[i]["label2"]  
        }) 

  return clusters1, clusters2
  
def assign_labels_to_all_clusters(node_dict, cluster, label):
  for i in xrange(len(node_dict)):
    for j in xrange(len(node_dict[i])): 
      if node_dict[i][j][label] == "none":
        node_dict[i][j][label] = node_dict[i][0][label]
        
  node_dict_list = make_clusters_into_single_list(node_dict)
  
  for i in xrange(len(cluster)):
    for j in xrange(len(node_dict_list)):
      if cluster[i]["id"] == node_dict_list[j]["id"]:
        return cluster 
