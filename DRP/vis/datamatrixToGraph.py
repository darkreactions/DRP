test_import json
from kdtree import *
from datamatrix import dataMatrix

# This will need to be modified to reflect current data (recommendations, etc.)
class graphNode:
  def __init__(self, name, friendlist, myID, headers, values, j):
    self.friendlist = friendlist
    self.inNodes = []
    self.dict = [values[headers[i]][j] for i in headers] 
    self.dict["name"] = name
    self.dict["id"] = myID 

  def getNumOut(self):
    return len(self.friendlist)

  def isEqual(self,otherNode):
    if (otherNode.dict["name"] == self.dict["name"] and otherNode.dict["id"] == self.dict["id"]):
      return True
    else:
      return False

import csv 
class myGraph:
  #Creates a graph from a datamatrix object as long as csv has: First col == id, Sec col == name of experiments; 3rd == inorg1
  # 4th col == inorg2, 5th col == org1, 6th(last) == purity; in between are fully numeric, filled in columns;
  def __init__(self, allreactions, headers, dids):
    
    headerValues = {}
    for i in xrange(len(headers)): 
        headersWithLists[headers[i]] = []

    edgeList = []    
    valuesArrs = []
   
    reader = csv.DictReader(datamatrix) 
    for row in reader:
        for i in xrange(len(headers)):
          headerValues[headers[i]].append(row[headers[i]])

    cleanedData = dataMatrix(headerValues) 
    cleanedData.removeCorrelatedLinregs(0,0)
          
    pointList = cleanedData.createPointList()
    tree = KDTree(pointList)
    for point in pointList:
      nnors = tree.findkNearestNeighbors(point, 5)
      #Because of the way index method works, will find only the first match (ignoring repeated elements)--hopefully not an issue  
      myFriendList = [pointList.index(neighbor) for neighbor in nnors] 
      edgeList.append(myFriendList)
    
    #datamatrix.dataset = np.transpose(datamatrix.dataset)
    #edgeList should be a list of lists (one list per point/row).

    #TODO: fix hardcoded parts 
    self.names = headervalues["names"]
    self.edgeList = edgeList
    self.ids = headerValues["reference"] 
    self.nodes = []
    self.numNodes = len(names)
    for i in range(0,self.numNodes):
        self.nodes.append(graphNode(self.names[i], self.edgeList[i], 1.0/self.numNodes, self.ids[i], headers, headerValues, i))

  def findNode(self, name):
    for node in self.nodes:
      if node.name == name:
        return node
    return -1

  def findNodebyID(self, testid):
    for node in self.nodes:
      if node.id == testid:
        return node
    return -1

  # formats node correctly for d3?
  def createDictForVis(self):
    visnodes = [] #TODO: This could be just self.nodes if we want ALL the data?
    for node in self.nodes:
      #all the headerValues need to be appended as a dictionary 
      visnodes.append({
        node.dict})

    dids = [] 
    for i in xrange(len(visnodes)):
     dids.append({ visnodes[i]["id"] }) 

    links = []
    for i in xrange(len(visnodes)):
      for target in self.edgeList[i]:
	    links.append({
            "source": i,
            "target": target,
            "value":1
          }
        )

    return {
        "nodes": visnodes,
        "links": links,
        "dids" : dids, 
  }
  def writeJson(self):
    import os
    from django.conf import settings
    self.setPageRanks()
    final_dict = self.createDictForVis()
    #f = open(os.path.join(os.path.dirname(os.path.realpath(__file__)), 'vis/matrix_formatted_for_vis'), "w")
    #dump = json.dump(final_dict, f, indent = 2)
    return final_dict  

  def sortNodesByPR(self):
    return sorted(self.nodes, key=lambda node: node.pagerank, reverse = True)


  def setInNodes(self):
    for node in self.nodes:
      node.inNodes = self.getInNodes(node)

  def getInNodes(self, node):
    inNodes = []
    for testnode in self.nodes:
      for testid in testnode.friendlist:
        if (testid == node.id):
          inNodes.append(testnode.id)
    return inNodes

  def numInNodes(self, node):
    return len(node.inNodes)

  def setPageRanks(self):
    self.setInNodes()
    counter = 0
    damping = 0.85
    converged = False #converges when sum of differences < .001
    while (not converged):
      ranklist = []
      for node in self.nodes:
        PR = (1-damping)/float(len(self.nodes))
        tsum = 0.0
        for friend in node.inNodes:
          currNode = self.findNodebyID(friend)
          if (currNode == -1):
            print "goddamnit..."
          if (self.numInNodes(currNode) == 0):
            tsum += currNode.pagerank/(len(self.nodes) - 1.0)
          else:
            tsum += (currNode.pagerank+0.0)/(currNode.getNumOut()+0.0)

        PR += tsum * damping
        ranklist.append(PR)

      rankdiff = 0
      for i in range(0,len(self.nodes)):
        diff = float(self.nodes[i].pagerank) - float(ranklist[i])
        rankdiff += diff
      self.updatePageRanks(ranklist)
      if (rankdiff < .000001):
        converged = True
      counter += 1

  def updatePageRanks(self, ranklist):
    for i in range(0,len(ranklist)):
      self.nodes[i].pagerank = ranklist[i]


  def getSinks(self):
    sinks = []
    for node in self.nodes:
      if len(node.friendlist) == 0:
        sinks.append(node)
    return sinks
