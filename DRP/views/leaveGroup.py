"""A view for a user to leave a group."""
from django.contrib.auth.decorators import login_required
from django.views.decorators.http import require_POST
from DRP.views.decorators import hasSignedLicense
from DRP.forms import LabGroupLeavingForm
from django.shortcuts import redirect

# TODO:change this so that memberships are tracked explicitly; at the moment this will break reaction validation if it is manually altered
# instead, change reaction validation so that it checks PRESENT members, and track membership as a boolean in an intermediary table
# between users and lab groups, this way validation will not break for
# historic members.


@login_required
@hasSignedLicense
@require_POST
def leaveGroup(request):
    """Allow a user to leave a group."""
    form = LabGroupLeavingForm(request.user, data=request.POST)
    if form.is_valid():
        request.user.labgroup_set.remove(form.cleaned_data.get('labGroup'))
    return redirect('account')
