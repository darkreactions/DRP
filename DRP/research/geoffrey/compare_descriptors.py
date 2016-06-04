#!/usr/bin/env python
"""Compare the values for two descriptors for the same reaction"""
import django
django.setup()
from DRP.models import BoolRxnDescriptor, OrdRxnDescriptor, NumRxnDescriptor, CatRxnDescriptor, BoolRxnDescriptorValue, OrdRxnDescriptorValue, NumRxnDescriptorValue, CatRxnDescriptorValue, PerformedReaction, DataSet, Compound
from sys import argv


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


def compare(d1_heading, d2_heading, comparison_function=None, desc1ValueType=BoolRxnDescriptorValue, desc2ValueType=BoolRxnDescriptorValue, reactions=None):
    if comparison_function is None:
        def _equal(v1, v2):
            return v1 == v2
        comparison_function = _equal
    if reactions is None:
        reactions = PerformedReaction.objects.all()

    d1_values = desc1ValueType.objects.filter(descriptor__heading=d1_heading)
    d2_values = desc2ValueType.objects.filter(descriptor__heading=d2_heading)
    print "found {} values for descriptor 1, {} for descriptor 2".format(d1_values.count(), d2_values.count())

    print d1_heading, d2_heading
    different = []

    count = 0
    for d1_val in d1_values:
        rxn = d1_val.reaction
        if reactions.filter(id=rxn.id).exists():
            count += 1
            d2_val = d2_values.get(reaction=rxn)
            if not comparison_function(d1_val.value, d2_val.value):
                if abs(d1_val.value - d2_val.value) != 1 or water not in rxn.compounds.all():
                    print rxn, d1_val.value, d2_val.value
                    print rxn.compounds.all()
                different.append(rxn)

    print "compared values for {} reactions".format(count)
    return different


def get_references(reactions):
    references = []
    for rxn in reactions:
        perf_rxn = PerformedReaction.objects.get(id=rxn.id)
        references.append(perf_rxn.reference)
    return references


def scale_function(scale, tol=0.1):
    def _scale(v1, v2):
        return (abs(scale * v1 - v2) < tol)
    return _scale


def shift_function(shift, tol=0.1):
    def _shift(v1, v2):
        return (abs(shift + v1 - v2) < tol)
    return _shift


def round_equal(v1, v2):
    return round(v1) == round(v2)


def bool_equal(v1, v2):
    return bool(v1) == bool(v2)

if __name__ == '__main__':
    scale60 = scale_function(60)
    shiftK = shift_function(273.15)
    tol_point_5 = shift_function(0, 0.55)
    #element = argv[1]
    water = Compound.objects.get(abbrev__iexact='H2O')

    d2_heading = 'boolean_valence_1'
    d1_heading = 'V0_legacy'
    reactions = DataSet.objects.get(name='valid_legacy_rxns_nonzero_compound').reactions.all()
    comparison_function = None
    desc1ValueType = getDescValueType(d1_heading)
    desc2ValueType = getDescValueType(d2_heading)
    different = compare(d1_heading, d2_heading, desc1ValueType=desc1ValueType, desc2ValueType=desc2ValueType, comparison_function=comparison_function, reactions=reactions)
    references = get_references(different)
    print len(references)
    print references
