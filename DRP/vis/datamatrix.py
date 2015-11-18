from get_data import *
from clustering import *
from kdtree import *

import linreg
import copy
import numpy as np
import statistics as stat
import math
import random
import sys
from kdtree import *

missing_data_token = "-1";

# Linreg = linear regression(actually correlation)
class dataMatrix:
  #Initialize data-matrix from dictionary of "header": [list of values corresponding to this header]
  # Properties of the datamatrix will be the header_list, dataset, num_cols, and num_rows
  def __init__(self,headersValuesDict):
    #calls get_data's method to find title of each column
    self.header_list = [key for key in headersValuesDict] 
    # sets number of columns equal to length of headers (length of the first row)--
    # initializes datamatrix
    self.dataset = headersValuesDict  # calls get_data method to put all rows of the file into a dictionary, with each cell in each row
    self.num_rows = len(self.dataset[header_list[0]]) #can pick any of the items in the header list
    self.num_cols = len(header_list)
    self.convertYesNotoOneZero()

  #Find the correlation between two columns of the data matrix
  def findLinreg(self, key1, key2):
    A = np.vstack([self.dataset[key1], np.ones(len(self.dataset[key2]))]).T
    m, b = np.linalg.lstsq(A, self.dataset[key2])[0]
    return m #, b

  #find all of the linregs. Returns a list of all comparisons in the form {Col1, Col2, Correl}
  def findAllLinregs(self):
    linreg_list = []
    for i in range(0,self.num_cols):
      for j in range(0, self.num_cols):
        if float(i) and float(j) and i != j:
          key1 = self.header_list[i]
          key2 = self.header_list[j]
          to_add = [self.header_list[i], self.header_list[j], self.findLinreg(key1, key2)]
          linreg_list.append(to_add)
    return linreg_list

  #Returns a list of all of the correlations in a list (useful for min/max linreg)
  def allLinregs(self):
    just_correls = []
    linreg_list = self.findAllLinregs()
    for item in linreg_list:
      just_correls.append(item[2])
    return just_correls 

'''
we don't want to do this
#Removes any linregs that are correlated within a certain range
  def removeCorrelatedLinregs(self, x, y):		
    correl_range = [.98,1.02]	
    for i in xrange(x, self.num_cols):
      for j in xrange(y, self.num_cols):
        if i != j and self.is_number(self.dataset[i]) and self.is_number(self.dataset[j]):
          if (-correl_range[0] > self.findLinreg(i,j) > -correl_range[1] or correl_range[0] < self.findLinreg(i,j) < correl_range[1]):
            if not (self.header_list[j] == "outcome"):
       	      self.removeColumn(j)
  	      return self.removeCorrelatedLinregs(i, j)
'''
  #Finds the mean value of a given column (helper for filling missing rows with mean)
  def meanCol(self,col_key):
    if (self.is_number(self.dataset[col_key][0]) == False):
      return 0
    total_sum = 0
    counter = 0
    for i in self.dataset[col_key]:
      num  = self.dataset[col_key][i]
      if num != missing_data_token and self.is_number(num): 
        total_sum += float(num)
      else:
        counter+=1
    return total_sum/(self.num_rows-counter)

  #Find standard deviation of a given column
  def stdevCol(self, col_key):
    if (self.is_number(self.dataset[col_key][0]) == False):
      return 0
    #find the variance
    mean = self.meanCol(col_key)
    squared_diffs = []
    #calculated all squared diffs
    for i in self.dataset[col_key]:
      num = self.dataset[col_key][i]
      try:
        squared_diffs.append((float(num)-mean)**2)
      except ValueError:
         pass 
	 #print "problem with"
        #print num
    #find the average of the squared_differences
    n = len(squared_diffs)
    sum_sdiffs = 0
    for item in squared_diffs:
      sum_sdiffs += item
    variance = sum_sdiffs/n
    return math.sqrt(variance)

  # Creates array of columns in matrix with normalized data
  def createStdevArrayAllCols(self):
    stdev_arr = [self.stdevCol(colkey) for colkey in self.headers]
    return stdev_arr

  #Creates array of means of all columns in data matrix (helpful when replacing missing tokens with mean)
  def createMeanArrayAllCols(self):
    mean_arr = [self.meanCol(colkey) for colkey in self.header_list]
    return mean_arr

  #Fills any missing data fields with the mean of its column. Missing data must be marked by a '?' or '-1'
  def fillMissingFieldsWithMean(self):
    counter = 0
    for i in xrange(0,self.num_cols):
      for j in xrange(0, self.num_rows):
        if self.dataset[i][j] == "?" or self.dataset[i][j] == missing_data_token:
          counter += 1
          self.dataset[i][j] = self.meanCol(i)
          self.cleaningflags[i][j] = 1

  #helper function used in removeStrings(), meanCol(), and stdevCol().
  def is_number(self, s):
    try:
        float(s)
        return True
    except:  
        return False

  #Convert yes's and no's to ones and zeroes
  def convertYesNotoOneZero(self):
    counter = 0
    for key in self.dataset:
        for j in xrange(len(self.dataset[key])): 
            if self.dataset[key][j] == "yes":
                counter += 1
                self.dataset[key][j] = 1
            elif self.dataset[key][j] == "no":
                counter += 1
                self.dataset[key][j] = 0
            elif self.dataset[key][j] == "?":
                counter += 1
                self.dataset[key][j] = 0

  #Creates array of mean values of each column, array of all normalized values (by column);
  #Creates matrix/array? of data points, mean of each column, and normalized data point
  def createPointList(self):
    meanArr = self.createMeanArrayAllCols()
    statArr = self.createStdevArrayAllCols()
    pointList = []
    for key in self.header_list:
            pointList.append(Point(self.dataset[key], meanArr, statArr))
    return pointList

  #FYI: Best not to use this function every time (slow)
  def findNN(self, QP):
    myTree = KDTree(self.createPointList())
    return myTree.findNearestNeighbor(QP)

