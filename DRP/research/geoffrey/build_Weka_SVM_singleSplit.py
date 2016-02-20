#!/usr/bin/env python

from DRP.models import PerformedReaction, ModelContainer, Descriptor
from django.db.models import Q
import operator
from sys import argv


def build_model(descriptor_header_file):
    reactions = PerformedReaction.objects.all()
    
    container = ModelContainer("weka", "SVM_PUK_basic", splitter="SingleSplitter",
                             reactions=reactions)
    container.save()
    
    
    #headers = ["reaction_temperature"]
    # get the headers to use from the descriptor header file
    with open(descriptor_header_file, 'r') as f:
        headers = [l.strip() for l in f.readlines()]
    
    predictors = get_descriptors_by_header(headers)
    responses = Descriptor.objects.filter(heading="boolean_crystallisation_outcome")
    
    container.build(predictors, responses)
    
    print container.summarize()

def get_descriptors_by_header(headers):
  Qs = [Q(heading=header) for header in headers]
  return Descriptor.objects.filter(reduce(operator.or_, Qs))

if __name__=='__main__':
  descriptor_header_file = argv[1]
  build_model(descriptor_header_file)
