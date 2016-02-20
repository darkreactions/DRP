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
    responses = Descriptor.objects.filter(heading="boolean_crystallisation_outcome")

    container.build(predictors, responses)

    conf_mtrcs = container.getConfusionMatrices()

    for model_mtrcs in conf_mtrcs:
        for descriptor_header, conf_mtrx in model_mtrcs:
            print "Confusion matrix for {}:".format(descriptor_header)
            print confusionMatrixString(conf_mtrx)


def confusionMatrixString(confusionMatrix, headers=True):
    """
    Returns a string that will display a confusionMatrix
    If headers=True, includes the headers as the first row and first column.
    """

    table = confusionMatrixTable(confusionMatrix, headers)
    return ('\n'.join([''.join(['{:^6}'.format(item) for item in row]) for row in table]))
    

def confusionMatrixTable(confusionMatrix, headers=True):
    """
    Converts a confusion matrix dictionary to a list of lists.
    Primarily for display purposes.
    Each list corresponds to a single true value and contains the
    counts for each predicted value.
    If headers=True, includes the headers as the first row and first column.
    """

    values = confusionMatrix.keys()
    table = [ [0 for predicted in values] for true in values] 

    for i, true in enumerate(values):
        for j, predicted in enumerate(values):
            table[j][i] = confusionMatrix[true][predicted]

    if headers:
        for j, predicted in enumerate(values):
            table[j].insert(0, str(predicted))
        table.insert(0, [""] + map(str, values))

    return table


def get_descriptors_by_header(headers):
    Qs = [Q(heading=header) for header in headers]
    return Descriptor.objects.filter(reduce(operator.or_, Qs))


if __name__ == '__main__':
    descriptor_header_file = argv[1]
    build_model(descriptor_header_file)
