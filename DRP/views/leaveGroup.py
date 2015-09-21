
from django.contrib.auth.decorators import login_required
from django.views.decorators.http import require_POST
from DRP.views.decorators import hasSignedLicense
from DRP.forms import LabGroupLeavingForm
from django.shortcuts import redirect

@login_required
@hasSignedLicense
@require_POST
def leaveGroup(request):
  form = LabGroupLeavingForm(request.user, data=request.POST)
  if form.is_valid():
    request.user.labgroup_set.remove(form.cleaned_data.get('labGroup'))
  return redirect('account')
  
