#!/usr/bin/env python

from DRP.models import PerformedReaction, ModelContainer, Descriptor, rxnDescriptorValues
from django.db.models import Q
import operator
import argparse

def build_model(reactions, predictors, responses, modelVisitorLibrary, modelVisitorTool, splitter):
    container = ModelContainer(modelVisitorLibrary, modelVisitorTool, splitter=splitter, reactions=reactions)
    container.save()
    container.build(predictors, responses)

    return container


def display_model_results(container):
    conf_mtrcs = container.getConfusionMatrices()

    for model_mtrcs in conf_mtrcs:
        for descriptor_header, conf_mtrx in model_mtrcs:
            print "Confusion matrix for {}:".format(descriptor_header)
            print confusionMatrixString(conf_mtrx)


def prepare_and_build_model(descriptor_headers, response_headers, modelVisitorLibrary, modelVisitorTool, splitter):
    #Grab all reactions with defined outcome descriptors
    reactions = PerformedReaction.objects.all()
    reactions = reactions.exclude(ordrxndescriptorvalue__in=rxnDescriptorValues.OrdRxnDescriptorValue.objects.filter(descriptor__heading__in=response_headers, value=None))
    reactions = reactions.exclude(boolrxndescriptorvalue__in=rxnDescriptorValues.BoolRxnDescriptorValue.objects.filter(descriptor__heading__in=response_headers, value=None))
    reactions = reactions.exclude(catrxndescriptorvalue__in=rxnDescriptorValues.CatRxnDescriptorValue.objects.filter(descriptor__heading__in=response_headers, value=None))

    predictors = Descriptor.objects.filter(heading__in=descriptor_headers)
    responses = Descriptor.objects.filter(heading__in=response_headers)

    container = build_model(reactions, predictors, responses, modelVisitorLibrary, modelVisitorTool, splitter)

    display_model_results(container)


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


#def get_descriptors_by_header(headers):
    #Qs = [Q(heading=header) for header in headers]
    #return Descriptor.objects.filter(reduce(operator.or_, Qs))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Builds a model', fromfile_prefix_chars='@')
    parser.add_argument('-p', '--predictor-headers', nargs='+')
    parser.add_argument('-r', '--response-headers', nargs='+', default=["boolean_crystallisation_outcome"])
    parser.add_argument('-ml', '--model-library', default="weka")
    parser.add_argument('-mt', '--model-tool', default="SVM_PUK_basic")
    parser.add_argument('-s', '--splitter', default="KFoldSplitter")
    args = parser.parse_args()
    
    ## get the headers to use from the descriptor header file
    #with open(args.descriptor_header_file, 'r') as f:
        #descriptor_headers = [l.strip() for l in f.readlines()]

    prepare_and_build_model(args.predictor_headers, args.response_headers, args.model_library, args.model_tool, args.splitter)
