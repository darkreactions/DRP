#!/usr/bin/env python
import django
django.setup()
from DRP.models import Descriptor
from DRP.models.rxnDescriptors import BoolRxnDescriptor, OrdRxnDescriptor, NumRxnDescriptor, CatRxnDescriptor
from DRP.models.predRxnDescriptors import PredBoolRxnDescriptor, PredOrdRxnDescriptor, PredNumRxnDescriptor, PredCatRxnDescriptor


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


def valid_descriptor_types():
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


def print_headers_with_names():
    for descriptor in valid_descriptors():
        print descriptor.heading, '\n\t', descriptor.name, '\n'

def print_headers():
    for descriptor in valid_descriptors():
        print descriptor.heading

def valid_descriptors():
    descriptor_list = []
    for descriptor in Descriptor.objects.all():
        descriptor_types = [BoolRxnDescriptor, OrdRxnDescriptor, NumRxnDescriptor, CatRxnDescriptor]
        predicted_descriptor_types = [PredBoolRxnDescriptor, PredOrdRxnDescriptor, PredNumRxnDescriptor, PredCatRxnDescriptor]
        
        if (is_descriptor_types(descriptor_types, descriptor) and
                not is_descriptor_types(predicted_descriptor_types, descriptor) and
                "outcome" not in descriptor.heading and 
                "rxnSpaceHash" not in descriptor.heading and 
                not descriptor.heading.startswith('_') and
                not descriptor.heading.startswith('transform') and 
                "examplepy" not in descriptor.heading):
            descriptor_list.append(descriptor)

    return descriptor_list

if __name__=='__main__':
    #print_headers()
    print_headers_with_names()
    #print_descriptor_types()
