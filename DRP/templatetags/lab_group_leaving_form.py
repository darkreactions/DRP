from django import template
from DRP.forms import LabGroupLeavingForm

register=template.Library()

def labGroupLeavingFormId():
  if not hasattr(compoundDeleteFormId, 'count'):
    compoundDeleteFormId.count = 0
  else:
    compoundDeleteFormId.count +=1
  return compoundDeleteFormId.count

@register.simple_tag(takes_context=True)
def lab_group_leaving_form(context, instance):
  return LabGroupLeavingForm(labGroup=instance, user=context['user'], auto_id='%s_{}'.format(compoundDeleteFormId())).as_ul()
