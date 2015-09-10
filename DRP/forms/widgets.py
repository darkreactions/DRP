'''Module containing some custom widgets'''

from django import forms
from django.utils import html

class SubmitButtonWidget(forms.Widget):
  '''Submit button widget borrowed from https://djangosnippets.org/snippets/2312'''

  def render(self, name, value, attrs="None"):
    return '<input type="submit" name="{}" value="{}" />'.format(html.escape(name), html.escape(value))
