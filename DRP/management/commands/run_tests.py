"""Command to run tests for DRP."""
from django.core.management.base import BaseCommand
from DRP.tests import suite, runTests
from django.conf import settings
import importlib
import unittest


class Command(BaseCommand):
    """Runs the tests using our custom testing framework."""

    help = 'Runs the full battery of DRP tests'

    def add_arguments(self, parser):
        """Add arguments for the argument parser."""
        parser.add_argument('test_modules', nargs='*',
                            help=('List of test modules to use. '
                                  'Should be dot-separated module paths starting at the test directory. '
                                  '(e.g. HttpTests.AboutPage to run just the tests for the about page. '
                                  'Default is all tests.'
                                  )
                            )
        parser.add_argument('--failfast', action='store_true',
                            help='Turn on the unittest failfast option. Tests halt at the first failure.'
                            )

    def handle(self, *args, **kwargs):
        """Handle the call for this command."""
        failfast = kwargs['failfast']
        if settings.TESTING:
            if kwargs['test_modules']:
                test_suite = unittest.TestSuite([getattr(importlib.import_module(
                    'DRP.tests.' + module), 'suite') for module in kwargs['test_modules']])
            else:
                test_suite = suite

            result = runTests(test_suite, failfast=failfast)
            if len(result.errors) > 0 or len(result.failures) > 0 or len(result.unexpectedSuccesses) > 0:
                exit(1)
        else:
            raise RuntimeError('Testing environment is not set')
