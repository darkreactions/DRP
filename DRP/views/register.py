'''A module containing only the register view'''
from DRP.forms import UserCreationForm
from django.template import RequestContext
from DRP.Email import Email
from django.shortcuts import render, redirect
from django.contrib.auth import login
from django.core.urlresolvers import reverse
from django.utils.http import urlencode

def register(request):
  '''A view to permit new users to sign up.'''
  redirected = False
  if request.method=='POST':
    form = UserCreationForm(request.POST)
    if form.is_valid():
      redirected = True
      form.save()
      if 'next' in request.GET.keys():
        return redirect('{0}?next={1}'.format(reverse('login'), urlencode(request.GET['next'])))
      else:
        return redirect('login')
  else:
    form = UserCreationForm()
  if not redirected:
    return render(request, 'register.html', RequestContext(request, {'form':form}))
