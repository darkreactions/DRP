"""
Miscellaneous utility functions for use in DRP
"""

from django.template.defaultfilters import slugify as _slugify

# When we move to python 3 we can use the builtin lrucache in functools
# There's a external library backport of functools from 3.2 to 2.7, but I opted to implement this function
# rather than adding another external dependency
def memoize(f):
    """
    Memoization decorator for a function taking one or more arguments.
    From here: http://code.activestate.com/recipes/578231-probably-the-fastest-memoization-decorator-in-the-/
    """
    class memodict(dict):
        def __getitem__(self, *key):
            return dict.__getitem__(self, key)

        def __missing__(self, key):
            self[key] = ret = f(*key)
            return ret

    return memodict().__getitem__

@memoize
def slugify(text):
    """Return a modified version of slug text.

    This modified version maintains compatibility with
    external languages such as R.
    """
    return _slugify(text).replace('-', '_')


@memoize
def generate_csvHeader(heading, calculatorSoftware, calculatorSoftwareVersion):
    """
    Method for generating csvHeader from Descriptor properties.
    This is a separate method to allow caching
    """
    return '{}_{}_{}'.format(
                             heading,
                             calculatorSoftware,
                             calculatorSoftwareVersion
                            )
    
