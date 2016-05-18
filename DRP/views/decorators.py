'''A group of decorators useful within the DRP views'''

from django.contrib.auth.decorators import login_required
from django.template import RequestContext
from django.http import HttpResponseForbidden
from django.template.loader import get_template
from DRP.models import LicenseAgreement, License
from DRP.models import LabGroup, PerformedReaction
from django.shortcuts import redirect
from django.core.urlresolvers import reverse
from django.utils.http import urlencode
from django.http import HttpResponseNotFound
from DRP.forms import LabGroupSelectionForm

def userHasLabGroup(view):
    '''This decorator checks that the user is a member of at least one lab group. Assumes login_required is an external decorator'''
    def hasLabGroup(request, *args, **kwargs):
        if not request.user.labgroup_set.all().exists():
            template = get_template('labgroup_403.html')
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
        if not License.objects.all().exists():
            template = get_template('license_404.html')
            return HttpResponseNotFound(template.render(RequestContext(request)))
        elif not LicenseAgreement.objects.filter(user=request.user, text=License.objects.latest()).exists():
            return redirect(reverse('license') + '?{0}'.format(urlencode({'next':request.path_info})))
        else:
            return view(request, *args, **kwargs)
    
    return _hasSignedLicense

def reactionExists(view, *args, **kwargs):
    """This decorator checks that a reaction exists before continuing with the internal view"""
    def _reactionExists(request, *args, **kwargs):
        rxn_id = kwargs['rxn_id']
        if PerformedReaction.objects.filter(id=rxn_id, labGroup__in=request.user.labgroup_set.all()).exists():
            return view(request, *args, **kwargs)
        else:
            raise Http404("This reaction cannot be found")

    return _reactionExists

def labGroupSelected(dispatch_method):
    '''Ensures a viewing lab group has been selected. This assumes a listview, hence it expects to decorate a method'''

    def _labGroupSelected(self, request, *args, **kwargs): 
        if request.user.is_authenticated():
            self.labForm = LabGroupSelectionForm(request.user)
            if request.user.labgroup_set.all().count() > 1:
                if 'labgroup_id' in request.session and request.user.labgroup_set.filter(pk=request.session['labgroup_id']).exists():
                    self.labGroup = request.user.labgroup_set.get(pk=request.session['labgroup_id'])
                    self.labForm.fields['labGroup'].initial = request.session['labgroup_id']
                elif 'labgroup_id' not in request.session:
                    return redirect(reverse('selectGroup') + '?{0}'.format(urlencode({'next':request.path_info})))
                else:
                    raise RuntimeError("This shouldn't happen")
                    self.labGroup = None
            elif request.user.labgroup_set.all().count() == 1:
                self.labGroup = request.user.labgroup_set.all()[0]
            else:
                self.labForm = None
                self.labGroup = None
        else:
            self.labForm = None
            self.labGroup = None
        return dispatch_method(self, request, *args, **kwargs)

    return _labGroupSelected
