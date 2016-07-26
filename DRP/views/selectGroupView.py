"""A module containing the view for selecting a viewing group on compound lists and reaction lists."""

from django.contrib.auth.decorators import login_required
from .decorators import hasSignedLicense, userHasLabGroup
from DRP.forms import LabGroupSelectionForm
from django.shortcuts import render, redirect
from django.template import RequestContext


@login_required
@hasSignedLicense
@userHasLabGroup
def selectGroup(request):
    """Make user select a group and then save that in the session. Deprecated and functionality to be removed."""
    if 'next' in request.GET:
        n = request.GET['next']
    else:
        n = '/database/'
    if request.method == "POST":
        form = LabGroupSelectionForm(request.user, data=request.POST)
        if form.is_valid():
            request.session['labgroup_id'] = form.cleaned_data['labGroup'].id
            return redirect(n)
        else:
            return render(request, 'select_group.html', {'form': form, 'next': n}, status=422)
    else:
        form = LabGroupSelectionForm(request.user)
        return render(request, 'select_group.html', {'form': form, 'next': n})
