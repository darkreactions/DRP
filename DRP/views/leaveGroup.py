
from django.views.decorators.http import require_POST
from DRP.forms import LabGroupLeavingForm
from django.shortcuts import redirect

@login_required
@hasSignedLicense
@require_POST
def leaveGroup(request)
  form = LabGroupLeavingForm(request.user)
  if form.is_valid():
    request.user.labgroup_set.remove(form.cleaned_data('labGroup'))
  return redirect('account')
  
