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

    d1_values = descValueType.objects.filter(descriptor__heading=d1_heading)
    d2_values = descValueType.objects.filter(descriptor__heading=d2_heading)
    print "found {} values for descriptor 1, {} for descriptor 2".format(d1_values.count(), d2_values.count())

    different = []
    if d1_values.count() < d2_values.count():
        for d1_val in d1_values:
            d2_val = d2_values.get(reaction=d1_val.reaction)
            if not comparison_function(d1_val.value, d2_val.value):
                print d1_val.reaction, d1_val.value, d2_val.value
                different.append(d1_val.reaction)
    else:
        for d2_val in d2_values:
            d1_val = d1_values.get(reaction=d2_val.reaction)
            if not comparison_function(d1_val.value, d2_val.value):
                print d1_val.reaction, d1_val.value, d2_val.value
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
        return (abs(scale*v1 - v2) < tol)
    return _scale

def shift_function(shift, tol=0.1):
    def _shift(v1, v2):
        return (abs(shift+v1 - v2) < tol)
    return _shift

def round_equal(v1, v2):
    return round(v1) == v2

if __name__ == '__main__':
    scale60 = scale_function(60)
    shiftK = shift_function(273.15)
    tol_point_5 = shift_function(0, 0.55)
    d1_heading = 'pH_legacy'
    d2_heading = 'reaction_pH'
    descValueType = getDescValueType(d1_heading)
    different = compare(d1_heading, d2_heading, descValueType=descValueType, comparison_function=round_equal)
    references = get_references(different)
    print len(references)
    print references
