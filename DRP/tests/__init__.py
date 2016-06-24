"""The suite of DRP tests."""
from . from drpTestCase import DRPTestCase, runTests
from . import fileTests
from . import unittest
from . import email
from . import forms
from . import httpTests
from . import decorators
from . import compoundDescriptor
from . import compoundToCsv
from . import compoundToArff
# import modelBuildingTests
# import DataImport
import modelValidators
# import splitters


suite = unittest.TestSuite([
    email.suite,
    forms.suite,
    httpTests.suite,
    # modelBuildingTests.suite,
    # DataImport.suite,
    compoundDescriptor.suite,
    compoundToCsv.suite,
    compoundToArff.suite,
    modelValidators.suite,
    # splitters.suite,
    fileTests.suite,
])


modules = [
    "email",
    "forms",
    "httpTests",
    # "modelBuildingTests",
    # "DataImport",
    "compoundFromCsv",
    "compoundDescriptor",
    "compoundToCsv",
    "compoundToArff",
    "fileTests",
    "modelValidators",
]
