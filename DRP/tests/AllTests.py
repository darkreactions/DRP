#!/usr/bin/python
import unittest, TestConfig
#As you create test files, add the import statements here
#import TemplateTest

suites = []
#suites = [TemplateTest.suite()]
#This list should be populated with tests 

if __name__ == '__main__':
    allTests = unittest.TestSuite(suites) 
    unittest.TextTestRunner(verbosity=2).run(allTests)

