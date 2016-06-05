'''Urls for pages related to authentication and user management'''
from django.conf.urls import url
import DRP.views
from django.contrib.auth.views import login, logout
from django.contrib.auth.decorators import login_required
from django.views.generic.base import TemplateView

urls = [
    url(r'^login.html$', login, {'template_name': 'login.html'}, name='login'),
    url(r'^logout.html$', logout, {'next_page': 'home'}),
    url(r'^register.html$', DRP.views.register, name='register'),
    url(r'^confirm.html$', DRP.views.confirm, name='confirm'),
    url(r'^license.html$', DRP.views.license, name='license'),
    url(r'^account/$', login_required(TemplateView.as_view(template_name='account.html')), name='account'),
    url(r'^account/leave_group$', DRP.views.leaveGroup, name='leaveGroup'),
    url(r'^account/join_group.html', DRP.views.joinGroup, name='joinGroup')
]
