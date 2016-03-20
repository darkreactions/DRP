#!/usr/bin/env python
import django
django.setup()
from DRP.models import Descriptor
from DRP.models import PredBoolRxnDescriptor, PredOrdRxnDescriptor, PredNumRxnDescriptor, PredCatRxnDescriptor
from DRP.models import PerformedReaction, DataSet, Descriptor, BoolRxnDescriptor, BoolRxnDescriptorValue, NumRxnDescriptor, NumRxnDescriptorValue, OrdRxnDescriptor, OrdRxnDescriptorValue, CatRxnDescriptor, CatRxnDescriptorValue
from django.db.models import Count
from sys import argv
import sys
from itertools import chain

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
    type_dict = {"Bool":0, "Ord":0, "Num":0, "Cat":0}
    for descriptor in valid_descriptors():
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
    
    for descriptor in descriptors:
        bqs = BoolRxnDescriptorValue.objects.filter(descriptor=descriptor, reaction__in=reactions)
        nqs = NumRxnDescriptorValue.objects.filter(descriptor=descriptor, reaction__in=reactions)
        oqs = OrdRxnDescriptorValue.objects.filter(descriptor=descriptor, reaction__in=reactions)
        cqs = CatRxnDescriptorValue.objects.filter(descriptor=descriptor, reaction__in=reactions)
        if bqs.exists():
            bqs = bqs.exclude(value=None)
            if bqs.exists():
                if len(bqs.values('value').annotate(num=Count('value'))) != 1: # the latter clause checks whether the values are unique
                    valid_descriptors.append(descriptor)
                else:
                    sys.stderr.write("{} excluded because all non-None values are the same\n".format(descriptor.heading))
            else:
                sys.stderr.write("{} excluded because all values are None\n".format(descriptor.heading))
        elif nqs.exists():
            nqs = nqs.exclude(value=None)
            if nqs.exists():
                if len(nqs.values('value').annotate(num=Count('value'))) != 1: # the latter clause checks whether the values are unique
                    valid_descriptors.append(descriptor)
                else:
                    sys.stderr.write("{} excluded because all non-None values are the same\n".format(descriptor.heading))
            else:
                sys.stderr.write("{} excluded because all values are None\n".format(descriptor.heading))
        elif oqs.exists():
            oqs = oqs.exclude(value=None)
            if oqs.exists():
                if len(oqs.values('value').annotate(num=Count('value'))) != 1: # the latter clause checks whether the values are unique
                    valid_descriptors.append(descriptor)
                else:
                    sys.stderr.write("{} excluded because all non-None values are the same\n".format(descriptor.heading))
            else:
                sys.stderr.write("{} excluded because all values are None\n".format(descriptor.heading))
        elif cqs.exists():
            cqs = cqs.exclude(value=None)
            if cqs.exists():
                if len(cqs.values('value').annotate(num=Count('value'))) != 1: # the latter clause checks whether the values are unique
                    valid_descriptors.append(descriptor)
                else:
                    sys.stderr.write("{} excluded because all non-None values are the same\n".format(descriptor.heading))
            else:
                sys.stderr.write("{} excluded because all values are None\n".format(descriptor.heading))
        else:
            sys.stderr.write("{} excluded because there are no values\n".format(descriptor.heading))
    return valid_descriptors

def print_headers_with_names(descriptors):
    for descriptor in descriptors:
        print descriptor.heading, '\n\t', descriptor.name, '\n'

def print_headers(descriptors):
    for descriptor in descriptors:
        print descriptor.heading

def valid_rxn_descriptors():
    descriptor_list = []
    for descriptor in chain(BoolRxnDescriptor.objects.all(), NumRxnDescriptor.objects.all(), 
                            OrdRxnDescriptor.objects.all(), CatRxnDescriptor.objects.all()):
        predicted_descriptor_types = [PredBoolRxnDescriptor, PredOrdRxnDescriptor, PredNumRxnDescriptor, PredCatRxnDescriptor]
        
        if (not is_descriptor_types(predicted_descriptor_types, descriptor) and
                "outcome" not in descriptor.heading and 
                "rxnSpaceHash" not in descriptor.heading and 
                not descriptor.heading.startswith('_') and
                not descriptor.heading.startswith('transform') and 
                "examplepy" not in descriptor.heading and
                "_legacy" not in descriptor.heading):
            descriptor_list.append(descriptor)

    return descriptor_list

if __name__=='__main__':
    #descs = valid_rxn_descriptors()
    ##print len(descs)
    #print_headers(descs)
    
    #print_headers()
    #print_headers_with_names()
    #print_descriptor_types()

     in_file = argv[1]
     reaction_set_name = argv[2]

     with open(in_file) as f:
         descriptor_headers = [l.strip() for l in f.readlines()]

     descriptors = Descriptor.objects.filter(heading__in=descriptor_headers)
     reactions = DataSet.objects.get(name=reaction_set_name).reactions.all()

     print_headers(filter_through_reactions(reactions, descriptors))
