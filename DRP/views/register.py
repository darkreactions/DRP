"""A module containing only the views for registration and confirmation."""
from DRP.forms import UserCreationForm
from django.template import RequestContext, Context
from DRP.Email import Email
from django.shortcuts import render, redirect
from django.contrib.auth import login
from django.core.urlresolvers import reverse
from django.template.loader import render_to_string
from DRP.forms import ConfirmationForm
from DRP.models import ConfirmationCode
from uuid import uuid4
from django.core.exceptions import PermissionDenied
from django.conf import settings


def register(request):
    """A view to permit new users to sign up."""

    if request.method == 'POST':
        form = UserCreationForm(request.POST)
        if form.is_valid():
            user = form.save()
            user.is_active = False
            user.save()
            code = uuid4()
            confCode = ConfirmationCode(user=user, code=code)
            confCode.save()
            emailText = render_to_string('confirmation_email.txt', {'name': user.first_name + ' ' + user.last_name, 'code': code, 'servername': settings.SERVER_NAME, 'testing': settings.TESTING})
            m = Email('Dark Reactions Project Registration', emailText, [user.email])
            m.send()
            return render(request, 'register.html', RequestContext(request, {'form': form, 'submitted': True}))
        else:
            return render(request, 'register.html', RequestContext(request, {'form': form, 'submitted': False}))
    else:
        form = UserCreationForm()
        return render(request, 'register.html', RequestContext(request, {'form': form, 'submitted': False}))


def confirm(request):
    """A view to confirm sign-up, and to render the license agreement binding."""

    if 'code' in request.GET.keys():
        if request.method == 'POST':
            form = ConfirmationForm(request, request.POST)
            if form.is_valid():
                user = form.get_user()
                if user.confirmationcode.code == request.GET['code']:
                    user.is_active = True
                    user.save()
                    return render(request, 'confirm.html', RequestContext(request, {'form': form, 'success': True, 'code': request.GET['code']}))
                else:
                    raise PermissionDenied()
            elif form.user_cache is not None:
                if form.user_cache.is_active:
                    return redirect(reverse('reactionlist'))
            else:
                return render(request, 'confirm.html', RequestContext(request, {'form': form, 'success': False, 'code': request.GET['code']}))
        else:
            form = ConfirmationForm(request)
            return render(request, 'confirm.html', RequestContext(request, {'form': form, 'success': False, 'code': request.GET['code']}))
    else:
        raise PermissionDenied()
