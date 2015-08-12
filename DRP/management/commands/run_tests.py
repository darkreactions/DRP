from django.core.management.base import BaseCommand
from DRP.tests import suite, runTests
from django.conf import settings

class Command(BaseCommand):
  help='Runs the full battery of DRP tests'''

  def handle(self, *args, **kwargs):
    if settings.TESTING:
      runTests(suite)
    else:
      raise RuntimeError('Testing environment is not set')
