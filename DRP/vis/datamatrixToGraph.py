

import json
from kdtree import *
from datamatrix import * 

# This will need to be modified to reflect current data (recommendations, etc.) 
class graphNode:
	def __init__(self, name, myID, friendlist, pagerank, purity, outcome, inorg1, inorg2, org1):
		self.name = name
		self.friendlist = friendlist
		self.id = int(myID)
		self.pagerank = pagerank
		self.purity = purity
		self.outcome = outcome
		self.inorg1 = inorg1
		self.inorg2 = inorg2
		self.org1 = org1
		#self.valuesArr = valuesArr
		self.inNodes = [] 
		
	def getNumOut(self):
		return len(self.friendlist)
		
	def isEqual(self,otherNode):
		if (otherNode.name == self.name and otherNode.id == self.id):
			return True
		else:
			return False

class myGraph:
	#Creates a graph from a datamatrix object as long as csv has: First col == id, Sec col == name of experiments; 3rd == inorg1
	# 4th col == inorg2, 5th col == org1, 6th(last) == purity; in between are fully numeric, filled in columns; 
	def __init__(self, datamatrix):
		idCol = 0
		namesCol = 1
		inorg1Col = 2
		inorg2Col = 3
		org1Col = 4
		linkslist = []	
		names = []
		ids = []
		valuesArrs = []
		inorg1s = []
		inorg2s = []
		org1s = []
		outcomes = []
		purities = []
		for r in range (datamatrix.num_rows):
			ids.append(datamatrix.dataset[idCol][r])
			names.append(datamatrix.dataset[namesCol][r])
			inorg1s.append(datamatrix.dataset[inorg1Col][r])
			inorg2s.append(datamatrix.dataset[inorg2Col][r])
			org1s.append(datamatrix.dataset[org1Col][r])
			outcomes.append(datamatrix.dataset[datamatrix.num_cols-1][r])
			purities.append(datamatrix.dataset[datamatrix.num_cols-2][r])
			
		
		datamatrix.removeColumn(1) #remove name
		datamatrix.removeColumn(1) #remove inorg1
		datamatrix.removeColumn(1) #remove inorg2
		datamatrix.removeColumn(1) #remove org1
		
		pointList = datamatrix.createPointList()
		tree = KDTree(pointList)
		
		for point in pointList:
			nnors = tree.findkNearestNeighbors(point, 5)
			myFriendList = []
			#print point.values[0]
			for neighbor in nnors:
				#print neighbor.values
				myFriendList.append(int(neighbor.values[0]))
			linkslist.append(myFriendList)
		
		datamatrix.removeColumn(0)		
		datamatrix.removeColumn(datamatrix.num_cols-1)
		datamatrix.removeColumn(datamatrix.num_cols-2)
		
		datamatrix.dataset = np.transpose(datamatrix.dataset)
		#to add dataset
		#for r in range (0, datamatrix.num_rows):
			#valuesArrs.append(datamatrix.dataset[r])
			
		self.names = names
		self.allLinks = linkslist
		self.nodes = []
		self.numNodes = len(names)
		#make graphNodes -- name, id, friendlist, pagerank, outcome)
		for i in range(0,self.numNodes):
			#print i
			self.nodes.append(graphNode(self.names[i], ids[i], self.allLinks[i], 1.0/self.numNodes, purities[i], outcomes[i], inorg1s[i], inorg2s[i], org1s[i]))
	
	
	
	def findNode(self, name):
		for node in self.nodes:
			#print node.name
			if node.name == name:
				return node
		return -1
	
	def findNodebyID(self, testid):
		for node in self.nodes:
			if node.id == testid:
				return node
		return -1
		
	def findNodePositionbyID(self, testid):
		for i in range (0, self.numNodes):
			if self.nodes[i].id == testid:
				return i
		return -1
	
	# formats node correctly for d3? 
	def createDictForVis(self):
		nodes = []
		for node in self.nodes:
			nodes.append({
			"name":node.name, 
			"pagerank":node.pagerank, 
			"id":node.id, 
			"purity":node.purity, 
			"outcome":node.outcome,
			"inorg1":node.inorg1,
			"inorg2":node.inorg2,
			"org1":node.org1
			}
			)
		links = []
		for source in range(0, self.numNodes):
			for target in self.allLinks[source]:
				links.append({
				"source": source, 
				"target": self.findNodePositionbyID(target), 
				"value":1
				}
				)
		return {
		"nodes": nodes, 
		"links": links}
	
	def writeJson(self):
		self.setPageRanks()
		final_dict = self.createDictForVis()
		out_file = raw_input('Write json file for FDG where? ...')
		with open(out_file, "wb") as fp:
			dump = json.dump(final_dict, fp, indent = 2)
			#out_file.close()
		return dump
	
	#returns list of nodes sorted by sortfield. Does not modify object.
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
			#print counter
			ranklist = []
			for node in self.nodes: 
				PR = (1-damping)/float(len(self.nodes))
				tsum = 0.0
				#print node.inNodes
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
			#print rankdiff
			#print ranklist
			self.updatePageRanks(ranklist)
			if (rankdiff < .000001):	
				converged = True
			counter += 1
			#print counter
			
	def updatePageRanks(self, ranklist):
		for i in range(0,len(ranklist)):
			self.nodes[i].pagerank = ranklist[i]
		
		
	def getSinks(self):
		sinks = []
		for node in self.nodes:
			if len(node.friendlist) == 0:
				sinks.append(node)
		return sinks			
	
	
		ru
# matrix = dataMatrix("filledNumericDataCleanwithTitleandId.csv")
# matrix.removeRandomRows(1000)
#print matrix.dataset[0]
#x = myGraph(matrix)
#x.setPageRanks()
#x.writeJson()

#print x.sortNodesByPR()[0].pagerank
#print x.sortNodesByPR()[200].pagerank

#print "A", x.nodes[0].pagerank
#print "B", x.nodes[1].pagerank
#print "C", x.nodes[2].pagerank
#print "D", x.nodes[3].pagerank

#print "names:", len(x.names)
#print "nodes:", x.nodes[1].pagerank



"""
def createNodes(input_file):
	nodelist = {"nodes":
	adjacencyList = open(input_file, 'r')
	for n in adjacencyList:
		nodelist += {"name": + n.name, "pagerank":n.pagerank}
	adjacencyList.close()

def createLinks(input_file):
	linklist = {"links":
	adjacencyList = open(input_file, 'r')
	for n in adjacencyList:
		for a in n:
			linklist += {"source": + n, "target": + a} 
	adjacencyList.close()

adjList = "myNetwork.json"
out_file = raw_input('Write json file where? ...')
outf = open(out_file, 'w')
final_dict = createNodes(adjList) + createLinks(adjList)
json.dump(graph_dict,outf)
outf.close()
"""
