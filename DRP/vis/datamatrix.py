
from get_data import *
from clustering import *
from kdtree import *

#import linreg
import copy
import numpy as np
#import statistics as stat
import math
import random
import sys
from kdtree import *

missing_data_token = "-1";

# Linreg = linear regression(actually crrelation)
#To loop through, go CR not RC (because rows and columns are flipped)
class dataMatrix:
  #Initialize data-matrix from file. Properties of the datamatrix will be the header_list,
  # dataset, num_cols, and num_rows
  def __init__(self,data):
    #calls get_data's method to find title of each column
    self.header_list = data.pop(0)
    # sets number of columns equal to length of headers (length of the first row)--
    # initializes datamatrix
    self.dataset = data  # calls get_data method to put all rows of the file into a dictionary, with each cell in each row
    self.num_rows = len(self.dataset)
    self.num_cols = len(self.dataset[0])
    self.convertYesNotoOneZero()

  #Automatically write matrix to a csv file, Filename is prese.t
  def writeToFile(self,filename):
    set_to_write = np.transpose(self.dataset)
    write_data(filename, self.header_list, set_to_write)

  # write a csv file of the normalized data
  def writeArray(self, array):
    out_file = raw_input('Write stdev with mean array where? ...')
    with open(out_file,'w') as f:
      f_csv = csv.writer(f)
      for i in range(len(array)):
        f_csv.writerow(array[i])

  # write file that, for every data point, indicates whether or not it has been modified(
  def writeCleaningFlagFile(self):
    out_file = raw_input('Write flag file where? ...')
    set_to_write = np.transpose(self.cleaningflags)
    write_data(out_file, self.header_list, set_to_write)

  #Find the correlation between two columns of the data matrix
  def findLinreg(self, col_num_one, col_num_two):
    A = np.vstack([self.dataset[col_num_one], np.ones(len(self.dataset[col_num_two]))]).T
    m, b = np.linalg.lstsq(A, self.dataset[col_num_two])[0]
    return m #, b

  #find all of the linregs. Returns a list of all comparisons in the form {Col1, Col2, Correl}
  def findAllLinregs(self):
    linreg_list = []
    for i in range(0,self.num_cols):
      for j in range(0, self.num_cols):
        if float(i) and float(j) and i != j:
          to_add = [self.header_list[i], self.header_list[j], self.findLinreg(i, j)]
          linreg_list.append(to_add)
    return linreg_list

  #Returns a list of all of the correlations in a list (useful for min/max linreg)
  def allLinregs(self):
    just_correls = []
    linreg_list = self.findAllLinregs()
    for item in linreg_list:
      just_correls.append(item[2])
    return just_

  #Finds the mean value of a given column (helper for filling missing rows with mean)
  def meanCol(self,num_col):
    if (self.is_number(self.dataset[0][num_col]) == False):
      return 0
    total_sum = 0
    counter = 0
    for row in self.dataset:
      i = row[num_col]
      if i != missing_data_token:
        total_sum += float(i)
      else:
        counter+=1
    return total_sum/(self.num_rows-counter)

  #Find standard deviation of a given column
  def stdevCol(self, num_of_col):
    if (self.is_number(self.dataset[0][num_of_col]) == False):
      return 0
    #find the variance
    mean = self.meanCol(num_of_col)
    squared_diffs = []
    #calculated all squared diffs
    for row in self.dataset:
      num = row[num_of_col]
      try:
        squared_diffs.append((float(num)-mean)**2)
      except ValueError:
        print "problem with"
        print num
    #find the average of the squared_differences
    n = len(squared_diffs)
    sum_sdiffs = 0
    for item in squared_diffs:
      sum_sdiffs += item
    variance = sum_sdiffs/n
    return math.sqrt(variance)

  # Creates array of columns in matrix with normalized data
  def createStdevArrayAllCols(self):
    stdev_arr = []
    for i in xrange(0, self.num_cols):
      stdev_arr.append(self.stdevCol(i))
    return stdev_arr

  #Creates array of means of all columns in data matrix (helpful when replacing missing tokens with mean)
  def createMeanArrayAllCols(self):
    mean_arr = []
    for i in xrange(0, self.num_cols):
      mean_arr.append(self.meanCol(i))
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
    except ValueError:
        return False

  #Convert yes's and no's to ones and zeroes
  def convertYesNotoOneZero(self):
    counter = 0
    for i in xrange(0,self.num_rows):
      for j in xrange(0, self.num_cols):
        if self.dataset[i][j] == "yes":
          counter += 1
          self.dataset[i][j] = 1
        elif self.dataset[i][j] == "no":
          counter += 1
          self.dataset[i][j] = 0
        elif self.dataset[i][j] == "?":
          counter += 1
          self.dataset[i][j] = 0

  #Remove non-numeric columns
  def removeStringCols(self):
    counter = 0
    for i in xrange(0,self.num_cols):
      for j in xrange(1, self.num_rows):
        if (self.is_number(self.dataset[i][j]) == False and not (self.dataset[i][j]) == "?"):
          counter += 1
          self.removeColumn(i)
          return self.removeStringCols()

  #??
  def isMadeUp(self, row, col):
    self.cleaningflags

  #Removes columns to speed up processing
  def removeRandomCols(self, numEndingCols):
    while (self.num_cols > numEndingCols):
      delCol = random.randint(10,self.num_cols-3)
      self.removeColumn(delCol)

  # Removes rows to speed up processing
  def removeRandomRows(self, numEndingRows):
    while (self.num_rows > numEndingRows):
      delRow = random.randint(0,self.num_rows-1)
      self.removeRow(delRow)

  #??
  def addIds(self):
    id_list = []
    for i in range(0, self.num_rows):
      id_list.append(i)
    col_list = [id_list]
    for c in range(0, self.num_cols):
      col_list.append(self.dataset[c])
    self.dataset = col_list
    new_headers = ["id"]
    for item in self.header_list:
      new_headers.append(item)
    self.header_list = new_headers

  #removes column
  def removeColumn(self, col_num):
    self.header_list.pop(col_num)
    self.dataset.pop(col_num)
    self.num_cols -= 1

  # Creates array of all headers and array of all columns
  def getOnlyCols(self, colarray):
    newDataset = []
    newHeaders = []
    for col in colarray:
      newHeaders.append(self.header_list[col])
      newDataset.append(self.dataset[col])

    #add purity and outcome
    newHeaders.append(self.header_list[self.num_cols-2])
    newDataset.append(self.dataset[self.num_cols-2])
    newHeaders.append(self.header_list[self.num_cols-1])
    newDataset.append(self.dataset[self.num_cols-1])
    self.dataset = newDataset
    self.header_list = newHeaders
    self.num_cols = len(self.header_list)
    self.num_rows = len(self.dataset[0])

  #Creates array of mean values of each column, array of all normalized values (by column);
  #Creates matrix/array? of data points, mean of each column, and normalized data point
  def createPointList(self):
    meanArr = self.createMeanArrayAllCols()
    statArr = self.createStdevArrayAllCols()
    pointList = []
    for row in self.dataset:
      pointList.append(Point(row, meanArr, statArr))
    return pointList

  #FYI: Best not to use this function every time (slow)
  def findNN(self, QP):
    myTree = KDTree(self.createPointList())
    return myTree.findNearestNeighbor(QP)

  #Prints the data matrix... poorly
  def __repr__(self):
    headers = ""
    data = ""
    for item in self.header_list:
      headers += item + ", "
    for item in self.dataset:
      data += str(item) + ", "
    return str(headers[:-1]) + "\n" + str(data[:-2])
