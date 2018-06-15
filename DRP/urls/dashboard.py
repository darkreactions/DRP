"""Dashboard's urls."""

from django.conf.urls import url
from django.views.generic.base import TemplateView
import DRP.views

urls = [
    url(r'^dashboard$', DRP.views.dashboard, name='dashboard')
]
