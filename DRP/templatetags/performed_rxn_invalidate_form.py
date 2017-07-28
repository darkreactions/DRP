"""Tags for generating Preformed reaction invalidation forms dynamically."""
from django import template
from DRP.forms import PerformedRxnInvalidateForm

register = template.Library()


def rxnInvalidateFormId():
    """Ensure that the html elements have unique id attrs."""
    if not hasattr(rxnInvalidateFormId, 'count'):
        rxnInvalidateFormId.count = 0
    else:
        rxnInvalidateFormId.count += 1
    return rxnInvalidateFormId.count


@register.simple_tag(takes_context=True)
def performed_rxn_invalidate_form(context, instance):
    """Generate the required form."""
    return PerformedRxnInvalidateForm(instance=instance, user=context['user'], auto_id='%s_invalidate_{}'.format(rxnInvalidateFormId())).as_ul()
