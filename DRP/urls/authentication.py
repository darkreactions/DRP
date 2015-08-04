'''Urls for pages related to authentication and user management'''
from django.conf.urls import patterns, include, url
    
urls = patterns('',
    url(r'^login.html$', login, {'template_name':'login.html'}, name='login'),
    (r'^logout.html$', logout, {'next_page':'home'}),
    url(r'^register.html$', DRP.views.register, name='register'),
    url(r'^confirm.html$', DRP.views.confirm, name='confirm'),
    url(r'^license.html$', DRP.views.license, name='license')
)
