'''This package provides tests model building in the DRP project'''
import unittest
import ModelContainer
import BCR
import J48
import KNN
import NaiveBayes

suite = unittest.TestSuite([
    ModelContainer.suite,
    BCR.suite,
    KNN.suite,
    J48.suite,
    NaiveBayes.suite,
])
