"""Compare the values for two descriptors for the same reaction"""
import django
django.setup()
from DRP.models import BoolRxnDescriptor, OrdRxnDescriptor, NumRxnDescriptor, CatRxnDescriptor, BoolRxnDescriptorValue, OrdRxnDescriptorValue, NumRxnDescriptorValue, CatRxnDescriptorValue, PerformedReaction


def getDescValueType(desc_heading):
    try:
        BoolRxnDescriptor.objects.get(heading=desc_heading)
        return BoolRxnDescriptorValue
    except BoolRxnDescriptor.DoesNotExist:
        pass
    try:
        OrdRxnDescriptor.objects.get(heading=desc_heading)
        return OrdRxnDescriptorValue
    except OrdRxnDescriptor.DoesNotExist:
        pass
    try:
        NumRxnDescriptor.objects.get(heading=desc_heading)
        return NumRxnDescriptorValue
    except NumRxnDescriptor.DoesNotExist:
        pass
    try:
        CatRxnDescriptor.objects.get(heading=desc_heading)
        return CatRxnDescriptorValue
    except CatRxnDescriptor.DoesNotExist:
        pass
    raise ValueError('Invalid heading given. Does not match any descriptor')


def compare(d1_heading, d2_heading, comparison_function=None, descValueType=BoolRxnDescriptorValue):
    if comparison_function is None:
        def _equal(v1, v2):
            return v1 == v2
        comparison_function = _equal 

    d1_values = BoolRxnDescriptorValue.objects.filter(descriptor__heading=d1_heading)
    d2_values = BoolRxnDescriptorValue.objects.filter(descriptor__heading=d2_heading)

    different = []
    if d1_values.count() < d2_values.count():
        for d1_val in d1_values:
            d2_val = d2_values.get(reaction=d1_val.reaction)
            if not comparison_function(d1_val.value, d2_val.value):
                different.append(d1_val.reaction)

    return different

def get_references(reactions):
    references = []
    for rxn in reactions:
        perf_rxn = PerformedReaction.objects.get(id=rxn.id)
        references.append(perf_rxn.reference)
    return references

def scale_function(scale, tol=0.1):
    def _scale(v1, v2):
        return (scale*v1 - v2 < tol)
    return scale_function

if __name__ == '__main__':
    scale60 = scale_function(60)
    d1_heading = 'boolean_outcome_legacy'
    d2_heading = 'boolean_crystallisation_outcome'
    descValueType = getDescValueType(d1_heading)
    different = compare(d1_heading, d2_heading, descValueType=descValueType)
    references = get_references(different)
    print references
