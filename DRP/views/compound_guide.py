from django.shortcuts import render
from django.contrib.auth.decorators import login_required

from DRP.models import *


#Send/receive the compound guide:
@login_required
def compound_guide(request):
  u = request.user
  lab_group = u.get_profile().lab_group
  guide = list(get_lab_CG(lab_group))
  print guide

  return render(request, 'compound_guide.html', {
    "guide": guide,
  })

