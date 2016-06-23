"""The suite of DRP tests."""
from DRPTestCase import DRPTestCase, runTests
import fileTests
import unittest
import Email
import forms
import HttpTests
import decorators
import CompoundDescriptor
import CompoundToCsv
import CompoundToArff
# import modelBuildingTests
# import DataImport
import modelValidators
# import splitters


suite = unittest.TestSuite([
    Email.suite,
    forms.suite,
    HttpTests.suite,
    # modelBuildingTests.suite,
    # DataImport.suite,
    CompoundDescriptor.suite,
    CompoundToCsv.suite,
    CompoundToArff.suite,
    modelValidators.suite,
    # splitters.suite,
    fileTests.suite,
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
