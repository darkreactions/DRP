'''A group of decorators useful within the DRP views'''

from django.contrib.auth.decorators import login_required
from django.template import RequestContext
from django.http import HttpResponseForbidden
from django.template.loader import get_template
from DRP.models import LicenseAgreement, License
from django.shortcuts import redirect
from django.core.urlresolvers import reverse
from django.utils.http import urlencode

def userHasLabGroup(view):
  '''This decorator checks that the user is a member of at least one lab group. Assumes login_required is an external decorator'''
  def hasLabGroup(request, *args, **kwargs):
    if not request.user.labgroup_set.all().exists():
      template=get_template('labgroup_403.html')
      return HttpResponseForbidden(template.render(RequestContext(request)))
    else:
      return view(request, *args, **kwargs)
  return hasLabGroup

def hasSignedLicense(view):
  '''This decorator checks that the user has signed the
  latest license and redirects them if not
  Assumes login_required is an external decorator
  '''
  def _hasSignedLicense(request, *args, **kwargs):
    if not LicenseAgreement.objects.filter(user=request.user, text=License.objects.latest()).exists():
      return redirect(reverse('license') + '?{0}'.format(urlencode({'next':request.path_info})))
    else:
      return view(request, *args, **kwargs)
  
  return _hasSignedLicense
