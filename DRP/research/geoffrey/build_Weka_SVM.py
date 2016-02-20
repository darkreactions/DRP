#!/usr/bin/env python

from DRP.models import PerformedReaction, ModelContainer, Descriptor
from django.db.models import Q
import operator
from sys import argv


def build_model(descriptor_header_file):
    reactions = PerformedReaction.objects.all()

    container = ModelContainer("weka", "SVM_PUK_basic", splitter="KFoldSplitter",
                                 reactions=reactions)
    container.save()

    #headers = ["reaction_temperature"]
    # get the headers to use from the descriptor header file
    with open(descriptor_header_file, 'r') as f:
        headers = [l.strip() for l in f.readlines()]

    predictors = get_descriptors_by_header(headers)
    responses = Descriptor.objects.filter(heading="crystallisation_outcome")

    container.build(predictors, responses)

    conf_mtrcs = container.getConfusionMatrices()
    #print len(conf_mtrcs)
    
    #for model_mtrcs in conf_mtrcs:
        #for descriptor_header, conf_mtrx in model_mtrcs:
            #print "Confusion matrix for {}:".format(descriptor_header)
            #print conf_mtrx


def display_confusion_matrix(confusion_matrix):
    """Display a confusion matrix dictionary as a table"""
    pass


def get_descriptors_by_header(headers):
    Qs = [Q(heading=header) for header in headers]
    return Descriptor.objects.filter(reduce(operator.or_, Qs))


if __name__ == '__main__':
    descriptor_header_file = argv[1]
    build_model(descriptor_header_file)
