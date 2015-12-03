'''Dashboard's urls'''

from django.conf.urls import patterns, include, url
from django.views.generic.base import TemplateView

urls = patterns('',
    url(r'^$', TemplateView.as_view(template_name="dashboard.html"), name='dashboard')
)
