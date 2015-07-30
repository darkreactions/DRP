#/usr/bin/env python
'''This module contains tests for the about page'''

import unittest
from HomePage import HomePage
loadTests = unittest.TestLoader().loadTestsFromTestCase

class AboutPage(HomePage):
  '''Performs a get request on the AboutPage to check for Html Validity'''

  url = HomePage.baseUrl + '/about.html'

  def setUp(self):
    self.response = requests.get(self.url)

suite = unittest.TestSuite([
  loadTests(AboutPage)
])

if __name__ == '__main__':
  #Runs the test- a good way to check that this particular test set works without having to run all the tests.
  unittest.main()
