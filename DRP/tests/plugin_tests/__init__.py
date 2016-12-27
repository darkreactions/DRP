"""The rdkit suite of tests for DRP."""
from . import drp_rdkit
import unittest

suite = unittest.TestSuite([
    drp_rdkit.suite
])

