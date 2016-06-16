"""A module of small useful functions for helpers."""
from django.shortcuts import redirect as django_redir
import urllib


def redirect(url, *args, **kwargs):
    """As per django's inbuilt redirect, but add get perameters to the uri."""
    params = kwargs.pop('params', {})
    response = django_redir(url, *args, **kwargs)
    if len(params) > 0:
        response['Location'] += '?' + urllib.urlencode(params)
    return response
