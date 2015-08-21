from django import template
from DRP.forms import CompoundDeleteForm

register = template.Library()

@register.simple_tag(takes_context=True)
def compound_delete_form(context, instance):
  return CompoundDeleteForm(instance=instance, user=context['user']).as_ul()
