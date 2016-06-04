
from django.shortcuts import redirect as django_redir
import urllib

def redirect(url, *args, **kwargs):
    params = kwargs.pop('params', {})
    response = django_redir(url, *args, **kwargs)
    if len(params) > 0:
        response['Location'] += '?' + urllib.urlencode(params)
    return response
