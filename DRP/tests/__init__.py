"""The suite of DRP tests."""
from .drpTestCase import DRPTestCase, runTests
from . import fileTests
import unittest
from . import email
from . import forms
from . import httpTests
from . import decorators
from . import compoundDescriptor
from . import compoundToCsv
from . import compoundToArff
# import modelBuildingTests
# import DataImport
from . import modelValidators
from . import plugin_tests
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
    plugin_tests.suite
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
