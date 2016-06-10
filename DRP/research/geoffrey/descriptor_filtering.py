#!/usr/bin/env python
import django
django.setup()
from DRP.models import DataSet, BoolRxnDescriptor, BoolRxnDescriptorValue, NumRxnDescriptor, NumRxnDescriptorValue, OrdRxnDescriptor, OrdRxnDescriptorValue, CatRxnDescriptor, CatRxnDescriptorValue
from DRP.models.querysets import MultiQuerySet
from django.db.models import Count
from sys import argv
import sys


def is_descriptor_type(descriptor_type, descriptor):
    """Returns boolean whether descriptor is of type descriptor_type"""
    try:
        descriptor_type.objects.get(id=descriptor.id)
        return True
    except:
        return False


def is_descriptor_types(descriptor_type_list, descriptor):
    """Returns boolean whether descriptor is any of the types in descriptor_type_list"""
    for descriptor_type in descriptor_type_list:
        if is_descriptor_type(descriptor_type, descriptor):
            return True
    return False


def descriptor_types(descriptors):
    type_dict = {"Bool": 0, "Ord": 0, "Num": 0, "Cat": 0}
    for descriptor in descriptors:
        if is_descriptor_type(BoolRxnDescriptor, descriptor):
            type_dict["Bool"] += 1
        if is_descriptor_type(OrdRxnDescriptor, descriptor):
            type_dict["Ord"] += 1
        if is_descriptor_type(CatRxnDescriptor, descriptor):
            type_dict["Cat"] += 1
        if is_descriptor_type(NumRxnDescriptor, descriptor):
            type_dict["Num"] += 1
    return type_dict


def print_descriptor_types():
    for descType, number in valid_descriptor_types().items():
        print "There are {} descriptors of the {} type".format(number, descType)


def filter_through_reactions(reactions, descriptors):
    valid_descriptors = []
    desc_val_types = [(NumRxnDescriptor, NumRxnDescriptorValue), (BoolRxnDescriptor, BoolRxnDescriptorValue),
                      (CatRxnDescriptor, CatRxnDescriptorValue), (OrdRxnDescriptor, OrdRxnDescriptorValue)
                      ]
    for descriptor in descriptors:
        for dt, vt in desc_val_types:
            if isinstance(descriptor, dt):
                qs = vt.objects.filter(descriptor=descriptor, reaction__in=reactions)

                if qs.exists():
                    qs = qs.exclude(value=None)
                    if qs.exists():
                        # the latter clause checks whether the values are unique
                        if qs.values('value').annotate(num=Count('value')).count() != 1:
                            valid_descriptors.append(descriptor)
                        else:
                            sys.stderr.write("{} excluded because all non-None values are the same\n".format(descriptor.heading))
                    else:
                        sys.stderr.write("{} excluded because all values are None\n".format(descriptor.heading))
                else:
                    sys.stderr.write("{} excluded because there are no values\n".format(descriptor.heading))
                break
    return valid_descriptors


def filter_through_reactions_nonMissing_unique(reactions, descriptors):
    desc_triples = []
    desc_val_types = [(BoolRxnDescriptor, BoolRxnDescriptorValue), (NumRxnDescriptor, NumRxnDescriptorValue),
                      (CatRxnDescriptor, CatRxnDescriptorValue), (OrdRxnDescriptor, OrdRxnDescriptorValue)
                      ]
    for descriptor in descriptors:
        for dt, vt in desc_val_types:
            if isinstance(descriptor, dt):
                qs = vt.objects.filter(descriptor=descriptor, reaction__in=reactions).exclude(value=None)

                nonNull = qs.count()
                distinct = qs.values('value').annotate(num=Count('value')).count()

                desc_triples.append((descriptor.heading, nonNull, distinct))

    return desc_triples


def print_headers_with_names(descriptors):
    for descriptor in descriptors:
        print descriptor.heading, '\n\t', descriptor.name, '\n'


def print_headers(descriptors, sort=False):
    headings = [d.heading for d in descriptors]
    if sort:
        headings.sort(key=lambda x: x.lower())
    print '\n'.join(headings)


def filter_qset(dqs, require_substring=[], require_prefix=[], require_suffix=[],
                    exclude_substring=[], exclude_prefix=[], exclude_suffix=[]):
    for s in exclude_substring:
        dqs = dqs.exclude(heading__contains=s)
    for s in exclude_prefix:
        dqs = dqs.exclude(heading__startswith=s)
    for s in exclude_suffix:
        dqs = dqs.exclude(heading__endswith=s)
    for s in require_substring:
        dqs = dqs.filter(heading__contains=s)
    for s in require_prefix:
        dqs = dqs.filter(heading__startswith=s)
    for s in require_suffix:
        dqs = dqs.filter(heading__endswith=s)

    return dqs


def rxn_descriptors():
    return MultiQuerySet(
                        NumRxnDescriptor.objects.all(),
                        BoolRxnDescriptor.objects.all(),
                        CatRxnDescriptor.objects.all(),
                        OrdRxnDescriptor.objects.all(),
                        )


def valid_rxn_descriptors():
    exclude_substring = ["_prediction_", "outcome", "rxnSpaceHash", "example"]
    exclude_prefix = ["_", "transform"]
    return filter_qset(rxn_descriptors(), exclude_substring=exclude_substring, exclude_prefix=exclude_prefix)


def valid_legacy_rxn_descriptors():
    require_suffix = ["_legacy", ]
    return filter_qset(valid_rxn_descriptors(), require_suffix=require_suffix)


def valid_nonlegacy_rxn_descriptors():
    exclude_suffix = ["_legacy", ]
    return filter_qset(valid_rxn_descriptors(), exclude_suffix=exclude_suffix)


def nonlegacy_pHless_rxn_descriptors():
    exclude_substring = ["_pH{}_".format(n) for n in range(0, 15)]
    return filter_qset(valid_nonlegacy_rxn_descriptors(), exclude_substring=exclude_substring)


def nonlegacy_nopHreaction_rxn_descriptors():
    exclude_substring = ["_pHreaction_",]
    return filter_qset(valid_nonlegacy_rxn_descriptors(), exclude_substring=exclude_substring)


def remove_CA(qset):
    exclude_substring = ["_chemaxoncxcalc_"]
    return filter_qset(qset, exclude_substring=exclude_substring)

if __name__ == '__main__':
    reaction_set_name = argv[1]

    descs = nonlegacy_nopHreaction_rxn_descriptors()
    descs = remove_CA(qsets)
    #descs = nonlegacy_pHless_rxn_descriptors()

    reactions = DataSet.objects.get(name=reaction_set_name).reactions.all()

    filtered_descriptors = filter_through_reactions(reactions, descs)
    sys.stderr.write("Kept {} of {} descriptors\n".format(len(filtered_descriptors), descs.count()))
    print_headers(filtered_descriptors, sort=True)
