'''A module containing the view for joining a group'''

from DRP.forms import LabGroupJoiningForm
from django.contrib.auth.decorators import login_required
from django.template.loader import get_template
from django.template import RequestContext
from DRP.models import LabGroup
from django.shortcuts import render
from django.http import HttpResponseNotFound

@login_required
def joinGroup(request):
  '''The view which governs the form for joining lab groups'''
  if LabGroup.objects.all().count() < 1:
    template = get_template('labgroup_404.html')
    return HttpResponseNotFound(template.render(RequestContext(request)))
  elif request.method=='POST':
    form = LabGroupJoiningForm(data=request.POST)
    if form.is_valid():
      form.cleaned_data['labGroup'].users.add(request.user)
      form = LabGroupJoiningForm()
    return render(request, 'join_group.html',{'form':form})
  else:
    form = LabGroupJoiningForm()
    return render(request, 'join_group.html', {'form':form})
