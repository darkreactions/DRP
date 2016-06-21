"""Context processors for the DRP project."""
from django.conf import settings


def testing(request):
    """
    Permit if blocks for filtering html that should only show when testing.

    This is used in particular for template testing.
    """
    if settings.TESTING:
        return {'testing': True}
    else:
        return {'testing': False}
