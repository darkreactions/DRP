
from get_data import *
#from clustering import *
#from datamatrix import *
import math


class Node:

    def __init__(self, point, dimension, parent):
        self.point = point
        self.left = None
        self.right = None
        self.parent = parent
        self.dimension = dimension

# for item in nnors:


class KDTree:

    def __init__(self, points):
        self.points = points
        self.root = self.makeKDTree(points, 0, None)

    def makeKDTree(self, pointList, splitIndex, parent):
        split_dimension = splitIndex % len(pointList)
        plen = len(pointList)
        if plen < 1:
            return None
        if plen == 1:
            return Node(pointList[0], split_dimension, parent)
        else:
            pointList = sorted(pointList, key=lambda point: point.values[split_dimension])
            median = int(plen / 2)
            node = Node(pointList[median], split_dimension, parent)
            node.left = self.makeKDTree(pointList[0:median], splitIndex + 1, node)
            node.right = self.makeKDTree(pointList[median:], splitIndex + 1, node)
            return node

    def findNearestNeighborHelper(self, search_point, current_node, current_best):
        if current_node.point.distTo(search_point) < current_best.point.distTo(search_point):
            current_best = current_node

        if current_node.left is None and current_node.right is None:
            return current_best
        elif current_node.right is None:
            current_best = self.findNearestNeighborHelper(search_point, current_node.left, current_best)
        elif (current_node.left is None):
            current_best = self.findNearestNeighborHelper(search_point, current_node.right, current_best)
        else:
            splitdim = current_node.dimension % current_node.point.d
            searchSplitVal = search_point.values[splitdim]
            searchVal = search_point.values[current_node.dimension]
            currentNodeSplitVal = current_node.point.values[splitdim]
            currentNodeVal = current_node.point.values[current_node.dimension]
            left = current_node.left
            right = current_node.right
            bestDist = search_point.distTo(current_best.point)
            try:
                if searchSplitVal < currentNodeSplitVal:
                    current_best = self.findNearestNeighborHelper(search_point, left, current_best)
                if (math.fabs(currentNodeVal - searchVal) < bestDist):
                    current_best = self.findNearestNeighborHelper(search_point, right, current_best)
                else:
                    current_best = self.findNearestNeighborHelper(search_point, right, current_best)
                    if math.fabs(currentNodeSplitVal - searchSplitVal) < bestDist:
                        current_best = self.findNearestNeighborHelper(search_point, left, current_best)
            except:
                left_best = self.findNearestNeighborHelper(search_point, left, current_best)
                right_best = self.findNearestNeighborHelper(search_point, right, current_best)
                leftDist = search_point.distTo(left_best.point)
                rightDist = search_point.distTo(right_best.point)
                tuples = [(leftDist, left_best), (rightDist, right_best), (bestDist, current_best)]
                current_best = min(tuples, key=lambda x: x[0])[1]
        return current_best

    def findNearestNeighbor(self, querypoint):
        return self.findNearestNeighborHelper(querypoint, self.root, self.root)

    def removePoint(self, point):
        newList = []
        for item in self.points:
            if (not item.equals(point)):
                newList.append(item)
        return newList

    def findkNearestNeighbors(self, querypoint, k):
        return self.kNearestNeighborsHelper(self.points, k, querypoint, [])

    def kNearestNeighborsHelper(self, points, k, QP, neighborList):
        if k == 0:
            return neighborList
        else:
            tree = KDTree(points)
            NN = tree.findNearestNeighbor(QP)
            # print NN.point.values[0:3]
            neighborList.append(NN.point)
            return self.kNearestNeighborsHelper(tree.removePoint(NN.point), k - 1, QP, neighborList)
