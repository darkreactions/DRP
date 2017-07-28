"""Tags for generating compound deletion forms dynamically."""
from django import template
from DRP.forms import CompoundDeleteForm

register = template.Library()


def compoundDeleteFormId():
    """Ensure that the html elements have unique id attrs."""
    if not hasattr(compoundDeleteFormId, 'count'):
        compoundDeleteFormId.count = 0
    else:
        compoundDeleteFormId.count += 1
    return compoundDeleteFormId.count


@register.simple_tag(takes_context=True)
def compound_delete_form(context, instance):
    """Generate a deletion form for a compound."""
    return CompoundDeleteForm(instance=instance, user=context['user'], auto_id='%s_delete_{}'.format(compoundDeleteFormId())).as_ul()
