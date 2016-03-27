"""
Miscellaneous utility functions for use in DRP
"""

from django.template.defaultfilters import slugify as _slugify
from functools32 import lru_cache #backport of pythion 3.2 functools with memoization decorator

def memodict(f):
    """ Memoization decorator for a function taking a single argument """
    class memodict(dict):
        def __missing__(self, key):
            ret = self[key] = f(key)
            return ret 
    return memodict().__getitem__


#@lru_cache(maxsize=64)
@memodict
def slugify(text):
    """Return a modified version of slug text.

    This modified version maintains compatibility with
    external languages such as R.
    """
    return _slugify(text).replace('-', '_')


#@lru_cache(maxsize=None)
def generate_csvHeader(heading, calculatorSoftware, calculatorSoftwareVersion):
    """
    Method for generating csvHeader from Descriptor properties.
    This is a separate method to allow caching
    """
    return '{}_{}_{}'.format(
                             heading,
                             slugify(calculatorSoftware),
                             calculatorSoftwareVersion
                            )
    
