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
    raise ValueError('Invalid heading {}. Does not match any descriptor'.format(desc_heading))


def compare(d1_heading, d2_heading, comparison_function=None, desc1ValueType=BoolRxnDescriptorValue, desc2ValueType=BoolRxnDescriptorValue, reactions=None):
    if comparison_function is None:
        def _equal(v1, v2):
            return v1 == v2
        comparison_function = _equal 
    if reactions is None:
        reactions = PerformedReaction.objects.all()

    d1_values = desc1ValueType.objects.filter(descriptor__heading=d1_heading, reaction__in=reactions)
    d2_values = desc2ValueType.objects.filter(descriptor__heading=d2_heading, reaction__in=reactions)
    print "found {} values for descriptor 1, {} for descriptor 2".format(d1_values.count(), d2_values.count())

    print d1_heading, d2_heading
    different = []

    if d1_values.count() > d2_values.count():
        for d1_val in d1_values:
            rxn = d1_val.reaction
            try:
                d2_val = d2_values.get(reaction=rxn)
            except BoolRxnDescriptorValue.DoesNotExist:
                print rxn.performedreaction, d1_val.value, 'DoesNotExist', rxn.compounds.all()
                different.append(rxn)
            else:
                if not comparison_function(d1_val.value, d2_val.value):
                    print rxn.performedreaction, d1_val.value, d2_val.value, rxn.compounds.all()
                    different.append(rxn)
    else:
        for d2_val in d2_values:
            rxn = d2_val.reaction
            try:
                d1_val = d1_values.get(reaction=rxn)
            except BoolRxnDescriptorValue.DoesNotExist:
                print rxn.performedreaction, d1_val.value, 'DoesNotExist', rxn.compounds.all()
                different.append(rxn)
            else:
                if not comparison_function(d1_val.value, d2_val.value):
                    print rxn.performedreaction, d1_val.value, d2_val.value, rxn.compounds.all()
                    different.append(rxn)

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
    return round(v1) == round(v2)

def bool_equal(v1, v2):
    return bool(v1) == bool(v2)

if __name__ == '__main__':
    scale60 = scale_function(60)
    shiftK = shift_function(273.15)
    tol_point_5 = shift_function(0, 0.55)
    #element = argv[1]
    water = Compound.objects.get(abbrev__iexact='H2O')
    
    for valence_num in range(0,8):
        d2_heading = 'Inorg_drpInorgAtom_boolean_valence_{}_DRP_1.5_any'.format(valence_num)
        d1_heading = 'V{}_legacy'.format(valence_num)
        reactions = DataSet.objects.get(name='valid_legacy_rxns_nonzero_compound').reactions.all()
        comparison_function = None
        desc1ValueType = getDescValueType(d1_heading)
        desc2ValueType = getDescValueType(d2_heading)
        different = compare(d1_heading, d2_heading, desc1ValueType=desc1ValueType, desc2ValueType=desc2ValueType, comparison_function=comparison_function, reactions=reactions)
        references = get_references(different)
        print len(references)
        print references
