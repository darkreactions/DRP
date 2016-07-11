"""A module containing the view for joining a group."""

from DRP.forms import LabGroupJoiningForm
from django.contrib.auth.decorators import login_required
from django.template.loader import get_template
from django.template import RequestContext
from DRP.models import LabGroup
from django.shortcuts import render
from django.http import HttpResponseNotFound
from .decorators import hasSignedLicense


@login_required
@hasSignedLicense
def joinGroup(request):
    """The view which governs the form for joining lab groups."""
    if not LabGroup.objects.all().exists():
        template = get_template('labgroup_404.html')
        return HttpResponseNotFound(template.render(RequestContext(request)))
    elif request.method == 'POST':
        form = LabGroupJoiningForm(data=request.POST)
        status = 200
        if form.is_valid():
            form.cleaned_data['labGroup'].users.add(request.user)
            form = LabGroupJoiningForm()
        else:
            status = 422
        return render(request, 'join_group.html', {'form': form}, status=422)
    else:
        form = LabGroupJoiningForm()
        return render(request, 'join_group.html', {'form': form})
