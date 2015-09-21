'''A module containing tests for the file structure and
standards conformance of the DRP'''

import unittest
import pep8
from django.conf import settings
import os.path
loadTests = unittest.TestLoader().loadTestsFromTestCase

_pep8Files = [
  ('DRP', 'settings_example.py')
  ('DRP', 'admin.py')
]

pep8Files = [os.path.join(*_pep8File) for _pep8File in _pep8Files]

class TestFiles(unittest.TestCase):
  '''The test class which tests the file structures'''

  def test_settingsExamplePresent(self):
    '''Test that the settings_example.py file is present'''
    self.assertTure(os.path.isfile(os.path.join(settings.BASE_DIR, 'settings_example.py'), 'settings example file is missing!')

  def testInitPresence(self):
    for fileName in pep8Files:
      fullFileName = os.path.join(settings.BASE_DIR, fileName)
      if os.path.isdir(fullFileName)
        self.assertTrue(os.path.isfile(os.path.join(fullFileName, '__init__.py')))

  def testPep8(self):
    for fileName in pep8Files:
      fullFileNames = [os.path.join(settings.BASE_DIR, fileName) for fileName in pep8Files]
      pep8style = pep8.StyleGuide()
      result = pep8Style.check_files(fullFileNames)
      self.assertEqual(results.total_errors, 0, 'Errors were found in pep8 conformance')

suite = unittest.TestSuite([
    loadTests(TestFiles)
  ])

if __name__=='__main__':
  unittest.TextTestRunner(verbosity=4).run(suite)
