"""Module containing some custom widgets."""

from django import forms
from django.utils import html


class SubmitButtonWidget(forms.Widget):
    """Submit button widget borrowed from https://djangosnippets.org/snippets/2312."""

    def __init__(self, value, attrs=None):
        """Gives the value (text) to this widget."""
        super(SubmitButtonWidget, self).__init__(attrs)
        self.value = value

    def render(self, name, value=None, attrs=None):
        """Does what you migth expect."""
        return '<input type="submit" name="{}" value="{}" />'.format(html.escape(name), html.escape(self.value))
