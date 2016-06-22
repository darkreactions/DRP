"""Tags for generating compound deletion forms dynamically."""
from django import template
from DRP.forms import PerformedRxnDeleteForm

register = template.Library()


def rxnDeleteFormId():
    """Ensure that the html elements have unique id attrs."""
    if not hasattr(rxnDeleteFormId, 'count'):
        rxnDeleteFormId.count = 0
    else:
        rxnDeleteFormId.count += 1
    return rxnDeleteFormId.count


@register.simple_tag(takes_context=True)
def performed_rxn_delete_form(context, instance):
    """Generate a deletion form for a performed reaction."""
    return PerformedRxnDeleteForm(instance=instance, user=context['user'], auto_id='%s_delete_{}'.format(rxnDeleteFormId())).as_ul()
