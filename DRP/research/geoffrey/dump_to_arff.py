#!/usr/bin/env python

from DRP.models import PerformedReaction, ModelContainer, Descriptor
from django.db.models import Q
import operator
from sys import argv


def dump_to_arff(descriptor_header_file, arff_file):
    reactions = PerformedReaction.objects.all()
    
    # get the headers to use from the descriptor header file
    with open(descriptor_header_file, 'r') as f:
        headers = [l.strip() for l in f.readlines()]
    
    descriptors = get_descriptors_by_header(headers)
    descriptor_headers = [d.csvHeader for d in descriptors]
    
    with open(arff_file, 'w') as f:
        reactions.toArff(f, expanded=True, whitelistHeaders=descriptor_headers)

def get_descriptors_by_header(headers):
    Qs = [Q(heading=header) for header in headers]
    return Descriptor.objects.filter(reduce(operator.or_, Qs))

if __name__=='__main__':
    descriptor_header_file = argv[1]
    arff_file = argv[2]
    dump_to_arff(descriptor_header_file, arff_file)
