

from math import *
from kdtree import *
import json
import random
from collections import Counter
# Returns the highest occurring item
#from datamatrix import *


class Point:

    def __init__(self, valuesArr, meanArr, statArr):
        self.id = valuesArr[0]
        self.values = valuesArr
        self.statArr = statArr
        self.meanArr = meanArr
        self.d = len(self.values)
        self.nv = self.normalized(self.values, self.meanArr, self.statArr)
        self.outcome = valuesArr[len(valuesArr) - 1]
        self.clusternum = 3000

    def normalized(self, valuesArr, meanArr, statArr):
        nv = []
        for i in range(0, len(valuesArr)):
            stdev = float(statArr[i])
            mean = float(meanArr[i])
            value = valuesArr[i]
            if stdev != 0 and mean != 0:
                nv.append(value - mean / stdev)
            else:
                nv.append(-9999)  # Just want to return a very large negative number
        return nv

    def equals(self, otherPoint):
        if (self.d != otherPoint.d):
            return False
        else:
            for i in range(0, self.d):
                if (self.values[i] != otherPoint.values[i]):
                    return False
        return True

    def distTo(self, otherPoint):
        mysum = 0
        for i in range(0, self.d):
            diff = self.nv[i] - otherPoint.nv[i]
            mysum += diff ** 2
        return sqrt(mysum)


class Cluster:

    def __init__(self, center):
        self.center = center
        self.points = []
        self.pointDistancesToCenter = []

    def addPoint(self, point, distToPoint):
        self.points.append(point)
        self.pointDistancesToCenter.append(distToPoint)

    def getDistToCenter(self, index):
        return self.pointDistancesToCenter[index]


class Clustering:

    def __init__(self, clusters):
        self.clusters = clusters
        self.k = len(clusters)

    def getCenters(self):
        centerList = []
        for i in range(0, self.k):
            centerList.append(self.clusters[i].center)
        return centerList

    def getAllPoints(self):
        pointList = []
        for c in range(0, self.k):
            for p in range(0, len(self.clusters[c].points)):
                self.clusters[c].points[p].clusternum = c
                # print "cnum GAP: ", self.clusters[c].points[p].clusternum
                pointList.append(self.clusters[c].points[p])
        return pointList


def findNodePosbyID(nodeList, testid):
    for i in range(0, len(nodeList)):
        if nodeList[i].id == testid:
            return i
    return -1


def findNodeID(nodeList, queryNode):
    for i in range(0, len(nodeList)):
        if nodeList[i].repPoint.equals(queryNode.repPoint):
            return nodeList[i].id
    return -1


def findNodePos(nodeList, queryNode):
    for i in range(0, len(nodeList)):
        if nodeList[i].repPoint.equals(queryNode.repPoint):
            return i
    return -1


def findPointbyID(pointList, testid):
    for i in range(0, len(pointList)):
        if pointList[i].id == testid:
            return pointList[i]
    return -1


class clusterNode:

    def __init__(self, repPoint, children):
        self.repPoint = repPoint
        self.id = self.repPoint.id
        self.children = children
        self.idsInCluster = [self.id]

    def addID(self, newID):
        self.idsInCluster.append(newID)

    def getPointsinCluster(self, allPoints):
        pointList = []
    # print self.idsInCluster
        for currID in self.idsInCluster:
            pointList.append(findPointbyID(allPoints, currID))
        return pointList

    def findNearestNode(self, clusterList):
        queryPoint = self.repPoint
        allPoints = []
        for node in clusterList:
            allPoints.append(node.repPoint)
        myTree = KDTree(allPoints)
        nearestPoint = myTree.findNearestNeighbor(queryPoint)
        nearestNode = clusterList[findNodePosbyID(clusterList, nearestPoint.point.id)]
        return nearestNode

    def equals(self, otherNode):
        return (self.id == otherNode.id and self.repPoint.equals(otherNode.repPoint))

    def isIn(self, lst):
        for item in lst:
            if self.equals(item):
                return True
        return False

    def findKNearestNodes(self, k, clusterList):
        print "FINDING K NEAREST NEIGHBORS!!!!!!!!!"
        return self.kNearestNeighborNodesHelper(k, clusterList, [], None)

    def kNearestNeighborNodesHelper(self, k, clusterList, neighborList, tree):
        if k == 0:
            threshhold = 2
            finalList = []
            bestDist = neighborList[0].repPoint.distTo(self.repPoint)

            for neighbor in neighborList:
                myDist = neighbor.repPoint.distTo(self.repPoint)
                if myDist < bestDist:
                    bestDist = myDist

            print "BestDist: ", bestDist
            for neighbor in neighborList:
                myDist = neighbor.repPoint.distTo(self.repPoint)

                if myDist > bestDist * threshhold or myDist == bestDist:
                    finalList.append(neighbor)
                    print "Dist: ", myDist, "Appended: True"
                else:
                    print "Dist: ", myDist, "Appended: False"

            return finalList
        else:
            if tree is None:
                allPoints = []
                for node in clusterList:
                    allPoints.append(node.repPoint)
                tree = KDTree(allPoints)
            QP = self.repPoint
            # print tree
            nearestPoint = tree.findNearestNeighbor(QP)

            # print nearestPoint.point.id
            # tree.removePoint(nearestPoint.point)
            nearestNode = clusterList[findNodePosbyID(clusterList, nearestPoint.point.id)]
            # print "nodeID: ", nearestNode.id

            neighborList.append(nearestNode)
            return self.kNearestNeighborNodesHelper(k - 1, clusterList, neighborList, KDTree(tree.removePoint(nearestPoint.point)))


def combineClusters(clusterOne, clusterTwo, allPoints):
    combinedCluster = clusterNode(clusterOne.repPoint, [clusterOne, clusterTwo])

    for item in clusterOne.idsInCluster:
        combinedCluster.idsInCluster.append(item)

    for item in clusterTwo.idsInCluster:
        combinedCluster.idsInCluster.append(item)

    combinedCluster.reLabel(allPoints)
    return combinedCluster


def combineClusterList(clusterList, allPoints):
    combineIndex = random.randint(0, len(clusterList) - 1)
    combinedCluster = clusterNode(clusterList[combineIndex].repPoint, clusterList)
    for cluster in clusterList:
        for item in cluster.idsInCluster:
            combinedCluster.idsInCluster.append(item)
    combinedCluster.reLabel(allPoints)
    return combinedCluster


def createHierCluster(points):
    clusterList = []
    for point in points:
        clusterList.append(clusterNode(point))

    while(len(clusterList) > 1):
        for i in range(0, len(clusterList) / 2):
            combineIndex = random.randint(0, len(clusterList) - 1)
            # print "CombInd ", combineIndex
            NN = clusterList[combineIndex].findNearestNode(clusterList)
            nodePos = findNodePos(clusterList, NN)

            # print "CLLen: ", len(clusterList), "CombInd: ", combineIndex, "NodePos: ", nodePos
            firstNode = clusterList.pop(combineIndex)
            if (nodePos >= combineIndex):
                nodePos = nodePos - 1
            secondNode = clusterList.pop(nodePos)
            clusterList.append(combineClusters(firstNode, secondNode, points))
    return clusterList[0]


def createHierClusterBetter(points):
    clusterList = []
    for point in points:
        clusterList.append(clusterNode(point))

    while (len(clusterList) > 1):
        # print len(clusterList)
        k = 20
        print "K ----- ", k
        if (len(clusterList) > k):
            combineIndex = random.randint(0, len(clusterList) - 1)
            NNors = clusterList[combineIndex].findKNearestNodes(k, clusterList)
            combineList = []
            for item in NNors:
                combineList.append(item)
            combineList.append(clusterList[combineIndex])
            clusterList.pop(combineIndex)
            for item in NNors:
                nodePos = findNodePos(clusterList, item)
                if (nodePos != -1):
                    print "CLLen: ", len(clusterList), "CombInd: ", combineIndex, "PopInd: ", nodePos, "id: ", item.id
                    clusterList.pop(nodePos)
            clusterList.append(combineClusterList(combineList, points))
        elif (len(clusterList) > 10):
            k = 10
            combineIndex = random.randint(0, len(clusterList) - 1)
            NNors = clusterList[combineIndex].findKNearestNodes(k, clusterList)
            combineList = []
            for item in NNors:
                combineList.append(item)
            combineList.append(clusterList[combineIndex])
            clusterList.pop(combineIndex)
            for item in NNors:
                nodePos = findNodePos(clusterList, item)
                if (nodePos != -1):
                    # print "CLLen: ", len(clusterList), "CombInd: ", combineIndex, "PopInd: ", nodePos, "id: ", item.id
                    clusterList.pop(nodePos)
            clusterList.append(combineClusterList(combineList, points))

        elif (len(clusterList) > 5):
            k = 5
            combineIndex = random.randint(0, len(clusterList) - 1)
            NNors = clusterList[combineIndex].findKNearestNodes(k, clusterList)
            combineList = []
            for item in NNors:
                combineList.append(item)
            combineList.append(clusterList[combineIndex])
            clusterList.pop(combineIndex)
            for item in NNors:
                nodePos = findNodePos(clusterList, item)
                if (nodePos != -1):
                    # print "CLLen: ", len(clusterList), "CombInd: ", combineIndex, "PopInd: ", nodePos, "id: ", item.id
                    clusterList.pop(nodePos)
            clusterList.append(combineClusterList(combineList, points))

        elif (len(clusterList) > 3):
            k = 3
            combineIndex = random.randint(0, len(clusterList) - 1)
            NNors = clusterList[combineIndex].findKNearestNodes(k, clusterList)
            combineList = []
            for item in NNors:
                combineList.append(item)
            combineList.append(clusterList[combineIndex])
            clusterList.pop(combineIndex)
            for item in NNors:
                nodePos = findNodePos(clusterList, item)
                if (nodePos != -1):
                    # print "CLLen: ", len(clusterList), "CombInd: ", combineIndex, "PopInd: ", nodePos, "id: ", item.id
                    clusterList.pop(nodePos)
            clusterList.append(combineClusterList(combineList, points))

        elif (len(clusterList) > 1):
            combineIndex = random.randint(0, len(clusterList) - 1)
            # print "CombInd ", combineIndex
            NN = clusterList[combineIndex].findNearestNode(clusterList)
            nodePos = findNodePos(clusterList, NN)

            firstNode = clusterList.pop(combineIndex)
            if (nodePos >= combineIndex):
                nodePos = nodePos - 1
            # print "CLLen: ", len(clusterList), "CombInd: ", combineIndex, "NodePos: ", nodePos
            if (nodePos != -1):
                secondNode = clusterList.pop(nodePos)
                clusterList.append(combineClusters(firstNode, secondNode, points))

    return clusterList[0]


def makeClustering(centers, points):
    clusters = []
    for c in range(0, len(centers)):
        newClust = Cluster(centers[c])  # could be a problem with PbyR?
        clusters.append(newClust)
    for p in range(0, len(points)):
        point = points[p]
        closestClusterIndex = 0
        closestClusterDist = centers[0].distTo(point)
        for c in range(0, len(centers)):
            dist = centers[c].distTo(point)
            if (dist < closestClusterDist):
                closestClusterDist = dist
                closestClusterIndex = c
            # print "cCI: ", closestClusterIndex
        clusters[closestClusterIndex].addPoint(point, closestClusterDist)
    return Clustering(clusters)


def writeHClusterJson(dictionary):
    final_dict = dictionary
    out_file = raw_input('Write json file for FDG where? ...')
    with open(out_file, "wb") as fp:
        dump = json.dump(final_dict, fp, indent=2)
    return dump
