"""Dashboard's urls."""

from django.conf.urls import url
from django.views.generic.base import TemplateView

urls = [
    url(r'^$', TemplateView.as_view(
        template_name="dashboard.html"), name='dashboard')
]
