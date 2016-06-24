"""The suite of DRP tests."""
from .drpTestCase import DRPTestCase, runTests
import .fileTests
import .unittest
import .email
import .forms
import .httpTests
import .decorators
import .compoundDescriptor
import .compoundToCsv
import .compoundToArff
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
