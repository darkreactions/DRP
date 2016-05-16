# !/usr/bin/env python
"""A module containing tests for standards conformance of the DRP."""

import unittest
import pep8
from django.conf import settings
import os
import sys
from cStringIO import StringIO
import pep257
loadTests = unittest.TestLoader().loadTestsFromTestCase

#_pep8Files = [
    #('DRP', 'settings_example.py'),
    ## directories and files seperated by commas,
    #('DRP', 'views', '__init__.py'),
    #('DRP', 'urls', 'public.py'),
    #('DRP', 'models', 'descriptors.py'),
    #('DRP', 'models', 'descriptorValues.py'),
    #('DRP', 'models', 'StatsModel.py'),
    #('DRP', 'models', 'Compound.py'),
    #('DRP', 'plugins', 'moldescriptors', 'example.py'),
    #('DRP', 'admin.py'),
    #('DRP', 'management', 'commands', 'check_hash_collisions.py'),
    ## you can finish on a directory test the whole thing,
    ## but this should be a python package if you do so
    ## (__init__.py must be present)
    ## ('DRP', 'admin.py'),
    #('DRP', 'tests', 'fileTests.py'),
    #('DRP', 'tests', 'CompoundDescriptor.py'),
    ## These commented out files are commented out because
    ## they would serve as examples had they been converted already!
#]

#pep8Files = [
    #os.path.join(
        #settings.BASE_DIR, os.path.join(*_pep8File)
    #) for _pep8File in _pep8Files
#]

#pep257Files = []

#for f in pep8Files:
    #if os.path.isdir(f):
        #for root, dirnames, fileNames in os.walk(f):
            #for fileName in fileNames:
                #if fileName.endswith('.py'):
                    #pep257Files.append(os.path.join(root, fileName))
    #else:
        #pep257Files.append(f)

settings_file = os.path.join(settings.APP_DIR, 'settings_example.py')


class OutputCapture:

    """A class as a workaround for broken parts in a pep8 module."""

    def __enter__(self):
        """Open the context manager, intercepting stdout."""
        self.text = ''
        self._stdout = sys.stdout
        sys.stdout = sys.sdout = self._stringIO = StringIO()
        return self

    def __exit__(self, *args, **kwargs):
        """Cease interception of sys.sdout."""
        self.text += self._stringIO.getvalue()
        sys.stdout = self._stdout

    def __str__(self):
        """String representation of the captured output."""
        return self.text


class TestFiles(unittest.TestCase):

    """Test that files conform to python standards."""

    def setUp(self):
        """Get list of all python files in project"""
        self.files = []
        for root, dirnames, fileNames in os.walk(settings.BASE_DIR):
            for fileName in fileNames:
                if fileName.endswith('.py'):
                    self.files.append(os.path.join(root, fileName))

    def test_pep8(self):
        """Test pep8 conformance."""
        with OutputCapture() as output:
            pep8style = pep8.StyleGuide(ignore=['E501'])
            result = pep8style.check_files(self.files)
        self.assertEqual(result.total_errors, 0,
                         "Found pep8 style errors (and warnings).\n{}".format(output))

    def test_pep257(self):
        """Test pep257 conformance."""
        # need to coerce into list to get the length.
        # Other option would be to write a generator has next function or something
        result = list(pep257.check(self.files))
        self.assertEqual(len(result), 0,
                         "Found pep257 errors (and warnings).\n{}".format('\n'.join([str(err) for err in result])))


    def test_settingsExamplePresent(self):
        """Test that the settings_example.py file is present."""
        self.assertTrue(
            os.path.isfile(settings_file),
            '{} is missing!'.format(settings_file)
        )

    def testInitPresence(self):
        """Test that any listed modules have the requisite __init__.py file."""
        for fileName in self.files:
            if os.path.isdir(fileName):
                self.assertTrue(
                    os.path.isfile(
                        os.path.join(fileName, '__init__.py')
                    )
                )

    #def testPep8(self):
        #"""Test that all listed files conform to pep8 standards."""
        #for fileName in pep8Files:
            #fullFileNames = [
                #os.path.join(
                    #settings.BASE_DIR, fileName
                #) for fileName in pep8Files
            #]
            #pep8Style = pep8.StyleGuide(ignore=['E501'])

            #with OutputCapture() as output:
                #result = pep8Style.check_files(fullFileNames)

                #self.assertEqual(
                    #result.total_errors,
                    #0,
                    #'{} errors were found in pep8 Conformance:\n{}'.format(result.total_errors, output)
                #)

    #def testPep257(self):
        #"""Test docstring conformance of all files."""
        ## errors = [str(e) for e in check(pep257Files)]
        ## self.assertEqual(len(errors), 0, '\n'.join(errors))
        #pass

suite = unittest.TestSuite([loadTests(TestFiles)])

if __name__ == '__main__':
    unittest.TextTestRunner(verbosity=4).run(suite)
