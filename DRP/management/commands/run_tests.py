from django.core.management.base import BaseCommand
from DRP.tests import suite, runTests
from django.conf import settings

class Command(BaseCommand):
  help='Runs the full battery of DRP tests'''

  def handle(self, *args, **kwargs):
    if settings.TESTING:
      result = runTests(suite)
      if len(result.errors) > 0 or len(result.failures)>0 or len(result.unexpectedSuccesses) > 0:
        exit(1)
    else:
      raise RuntimeError('Testing environment is not set')
