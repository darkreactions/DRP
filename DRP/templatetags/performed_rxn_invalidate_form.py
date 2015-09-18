
from django import template
from DRP.forms import PerformedRxnInvalidateForm

register = template.Library()

def rxnInvalidateFormId():
  if not hasattr(rxnDeleteFormId, 'count'):
    rxnDeleteFormId.count = 0
  else:
    rxnDeleteFormId.count +=1
  return rxnDeleteFormId.count

@register.simple_tag(takes_context=True)
def performed_rxn_invalidate_form(context, instance):
  return PerformedRxnInvalidateForm(instance=instance, user=context['user'], auto_id='%s_{}'.format(rxnInvalidateFormId())).as_ul()
