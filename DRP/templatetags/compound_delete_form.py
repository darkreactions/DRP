from django import template
from DRP.forms import CompoundDeleteForm

register = template.Library()

@register.simple_tag
def compound_delete_form(instance):
  return CompoundDeleteForm(instance=instance).as_ul()
