#!/usr/bin/env python

from sys import argv
import operator
from DRP.research.geoffrey.distance_learning.metricLearn import LMNN
from DRP.models import PerformedReaction, ModelContainer, Descriptor, rxnDescriptorValues
from django.db.models import Q
from functools import reduce


def train(descriptor_header_file, outfile):
    reactions = PerformedReaction.objects.all()

    # get the headers to use from the descriptor header file
    with open(descriptor_header_file, 'r') as f:
        descriptorHeaders = [l.strip() for l in f.readlines()]

    predictors = get_descriptors_by_header(descriptorHeaders)
    responses = Descriptor.objects.filter(
        heading="boolean_crystallisation_outcome")

    predictorHeaders = [d.csvHeader for d in predictors]
    responseHeaders = [d.csvHeader for d in responses]

    lmnn = LMNN(reactions, predictorHeaders, responseHeaders)

    lmnn.train()

    with open(outfile, 'wb') as f:
        lmnn.save(f)


def get_descriptors_by_header(headers):
    Qs = [Q(heading=header) for header in headers]
    return Descriptor.objects.filter(reduce(operator.or_, Qs))

if __name__ == '__main__':
    descriptor_header_file = argv[1]
    outfile = argv[2]
    train(descriptor_header_file, outfile)
