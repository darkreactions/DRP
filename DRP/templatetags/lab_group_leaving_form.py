from django import template
from DRP.forms import LabGroupLeavingForm

register = template.Library()


def labGroupLeavingFormId():
    if not hasattr(labGroupLeavingFormId, 'count'):
        labGroupLeavingFormId.count = 0
    else:
        labGroupLeavingFormId.count += 1
    return labGroupLeavingFormId.count


@register.simple_tag(takes_context=True)
def lab_group_leaving_form(context, instance):
    return LabGroupLeavingForm(labGroup=instance, user=context['user'], auto_id='%s_{}'.format(labGroupLeavingFormId())).as_ul()
