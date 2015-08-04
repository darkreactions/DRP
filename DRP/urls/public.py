'''Publicly viewable pages'''

from django.conf.urls import patterns, include, url

urls = patterns('',
    url(r'^$', TemplateView.as_view(template_name="home.html"), name='home'),
    (r'^about.html$', TemplateView.as_view(template_name="about.html")),
    url(r'^contact.html$', DRP.views.contact, name='contact')
)

