'''This package provides tests model building in the DRP project'''
import unittest
import J48
import KNN
import NaiveBayes
import SVM_basic
import SVM_BCR

suite = unittest.TestSuite([
    J48.suite,
    KNN.suite,
    NaiveBayes.suite,
    SVM_basic.suite,
    SVM_BCR.suite,
])
