'''A group of decorators useful within the DRP views'''

from django.contrib.auth.decorators import login_required
from django.template import RequestContext
from django.http import HttpResponseForbidden
from django.template.loader import get_template

def userHasLabGroup(view):
  '''This decorator checks that the user is logged in
  AND is a member of at least one lab group
  '''
  def hasLabGroup(request, *args, **kwargs):
    if request.user.labgroup_set.all().count() < 1:
      template=get_template('labgroup_403.html')
      return HttpResponseForbidden(template.render(request, RequestContext(request)))      
    else:
      return view(request, *args, **kwargs)
  return login_required(hasLabGroup)
