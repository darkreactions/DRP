#!/usr/bin/env python

from DRP.models import PerformedReaction, ModelContainer, Descriptor, rxnDescriptorValues
from sys import argv
from DRP.research.geoffrey.distance_learning.metricLearn import ITML

def train(descriptor_header_file):
    reactions = PerformedReaction.objects.all()[:100]
  
    # get the headers to use from the descriptor header file
    with open(descriptor_header_file, 'r') as f:
        descriptorHeaders = [l.strip() for l in f.readlines()]

    responseHeaders = ["boolean_crystallisation_outcome"]
    itml = ITML()

    itml.train(reactions, descriptorHeaders, responseHeaders)


if __name__=='__main__':
  descriptor_header_file = argv[1]
  train(descriptor_header_file)
