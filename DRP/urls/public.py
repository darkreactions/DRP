"""Publicly viewable pages."""

from django.conf.urls import url
from django.views.generic.base import TemplateView
import DRP.views

urls = [
    url(r'^$', TemplateView.as_view(template_name="home.html"), name='home'),
    url(r'^about.html$', TemplateView.as_view(template_name="about.html")),
    url(r'^contact.html$', DRP.views.contact, name='contact'),
]
