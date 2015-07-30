'''A module containing only the views for registration and confirmation'''
from DRP.forms import UserCreationForm
from django.template import RequestContext
from DRP.Email import Email
from django.shortcuts import render, redirect
from django.contrib.auth import login
from django.core.urlresolvers import reverse
from django.utils.http import urlencode
from django.template.loader import render_to_string
from drp.forms import ConfirmationForm
from DRP.models import ConfirmationCode
from uuid import uuid4
from django.core.exceptions import PermissionDenied

def register(request):
  '''A view to permit new users to sign up.'''
  
  if request.method=='POST':
    form = UserCreationForm(request.POST)
    if form.is_valid():
      user = form.save()
      user.active = False
      user.save()
      code = uuid4()
      confCode = ConfirmationCode(user=user, code=code)
      confCode.save()
      emailText = render_to_string('confirmation_email.html', {'name':user.first_name + ' ' + user.second_name, 'code':code})
      Email('Dark Reactions Project Registration', emailText, [user.email])
      return render(request, 'register.html', {'form':form, 'submitted':True})
    else:
      return render(request, 'register.html', {'form':form, 'submitted':False})
  else:
    form = UserCreationForm()
    return render(request, 'register.html', RequestContext(request, {'form':form, 'submitted':False}))

def confirm(request):
  '''A view to confirm sign-up, and to render the license agreement binding'''

  if 'code' in request.GET.keys():
    if request.method = 'POST':
      form = ConfirmationForm(request, request.POST) 
      if form.is_valid():
        user = form.get_user()
        if user.confirmationcode = request.GET['code']:
          user.active = True
          user.save()
          return render(request, 'confirm.html', RequestContext(request, {'form':form, 'success':True, 'code':urlencode(request.GET['code']})))
        else:
          raise PermissionDenied()
      else:
        return render(request, 'confirm.html', RequestContext(request, {'form':form, 'success':False, 'code':urlencode(request.GET['code']}))
    else:
      form = ConfirmationForm(request)
      render(request, 'confirm.html', RequestContext(request, {'form':form, 'success':False, 'code':urlencode(request.GET['code']}))
  else:
    raise PermissionDenied()
