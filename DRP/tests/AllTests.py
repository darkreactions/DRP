#!/usr/bin/env python
import unittest
#As you create test files, add the import statements here
#import TemplateTest

suites = []
#suites = [TemplateTest.suite()]
#This list should be populated with tests 

if __name__ == '__main__':
    allTests = unittest.TestSuite(suites) 
    results = unittest.TestResult()
    allTests.run(results)
    if len(results.errors == 0) and len(results.failures) == 0 and len(results.unexpectedSuccesses) == 0:
        print('All tests passed!')
        exitStatus = 0
    else:
        print('Not all Tests Passed. Generating Report:\n')
        exitStatus = 1
        if len(results.errors) > 0:
            print('The following tests raised unexpected exceptions:\n')
            for e in results.errors:
                print(e[0].id())
        if len(results.failures) > 0:
            for f in results.failures:
                print(f[0].id())
        if len(results.unexpectedSuccesses) > 0
            for u in results.unexpectedSuccesses:
                print(u.id())
    if len(results.skipped) > 0:
        print('The following tests were skipped:')
        for s in results.skipped:
            print(s[0].id())
    exit(exitStatus)
