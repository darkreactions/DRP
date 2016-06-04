from DRPTestCase import DRPTestCase, runTests
import fileTests
import unittest
import Email
import forms
import HttpTests
import decorators
import CompoundFromCsv
import CompoundDescriptor
import CompoundToCsv
import CompoundToArff
#import modelBuildingTests
#import DataImport
import modelValidators


suite = unittest.TestSuite([
    Email.suite,
    forms.suite,
    HttpTests.suite,
    #    modelBuildingTests.suite,
    #    DataImport.suite,
    CompoundFromCsv.suite,
    CompoundDescriptor.suite,
    CompoundToCsv.suite,
    CompoundToArff.suite,
    fileTests.suite,
    modelValidators.suite,
])


modules = [
    "Email",
    "forms",
    "HttpTests",
    # "modelBuildingTests",
    # "DataImport",
    "CompoundFromCsv",
    "CompoundDescriptor",
    "CompoundToCsv",
    "CompoundToArff",
    "fileTests",
    "modelValidators",
]
