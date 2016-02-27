#!/usr/bin/env python

import django
django.setup()
from DRP.models import PerformedReaction, ModelContainer, Descriptor, rxnDescriptorValues
import operator
import argparse


def build_model(reactions, predictors, responses, modelVisitorLibrary, modelVisitorTool, splitter, verbose=False):
    container = ModelContainer.create(modelVisitorLibrary=modelVisitorLibrary, modelVisitorTool=modelVisitorTool, splitter=splitter, reactions=reactions)
    container.save()
    container.build(predictors, responses, verbose=verbose)
    container.save()

    return container


def display_model_results(container):
    conf_mtrcs = container.getConfusionMatrices()

    for model_mtrcs in conf_mtrcs:
        for descriptor_header, conf_mtrx in model_mtrcs:
            print "Confusion matrix for {}:".format(descriptor_header)
            print confusionMatrixString(conf_mtrx)
            print "Accuracy: {:.3}".format(accuracy(conf_mtrx))
            print "BCR: {:.3}".format(BCR(conf_mtrx))


def prepare_build_display_model(descriptor_headers, response_headers, modelVisitorLibrary, modelVisitorTool, splitter, verbose=False):
    """
    Build and display a model with the specified tools
    """
    # Grab all reactions with defined outcome descriptors
    reactions = PerformedReaction.objects.all()
    reactions = reactions.exclude(ordrxndescriptorvalue__in=rxnDescriptorValues.OrdRxnDescriptorValue.objects.filter(descriptor__heading__in=response_headers, value=None))
    reactions = reactions.exclude(boolrxndescriptorvalue__in=rxnDescriptorValues.BoolRxnDescriptorValue.objects.filter(descriptor__heading__in=response_headers, value=None))
    reactions = reactions.exclude(catrxndescriptorvalue__in=rxnDescriptorValues.CatRxnDescriptorValue.objects.filter(descriptor__heading__in=response_headers, value=None))

    predictors = Descriptor.objects.filter(heading__in=descriptor_headers)
    responses = Descriptor.objects.filter(heading__in=response_headers)

    container = build_model(reactions, predictors, responses, modelVisitorLibrary, modelVisitorTool, splitter, verbose=verbose)

    display_model_results(container)

def accuracy(conf):
    correct = 0.0
    total = 0.0
    for true, guesses in conf.items():
        for guess, count in guesses.items():
            if true == guess:
                correct += count
            total += count
    return correct/total

    
def BCR(conf):
    class_accuracy_sum = 0.0
    num_classes = 0.0
    for true, guesses in conf.items():
        class_correct = 0.0
        class_total = 0.0
        for guess, count in guesses.items():
            if true == guess:
                class_correct += count
            class_total += count
        if class_total != 0:
            class_accuracy_sum += class_correct/class_total
            num_classes += 1
    
    return class_accuracy_sum/num_classes

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
    table = [[0 for predicted in values] for true in values]

    for i, true in enumerate(values):
        for j, predicted in enumerate(values):
            table[j][i] = confusionMatrix[true][predicted]

    if headers:
        for j, predicted in enumerate(values):
            table[j].insert(0, str(predicted))
        table.insert(0, [""] + map(str, values))

    return table

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Builds a model', fromfile_prefix_chars='@',
                                     epilog="Prefix arguments with '@' to specify a file containing newline"
                                     "-separated values for that argument. e.g.'-p @predictor_headers.txt'"
                                     " to pass multiple descriptors from a file as predictors")
    parser.add_argument('-p', '--predictor-headers', nargs='+',
                        help='One or more descriptors to use as predictors.', required=True)
    parser.add_argument('-r', '--response-headers', nargs='+', default=["boolean_crystallisation_outcome"],
                        help='One or more descriptors to predict. '
                        'Note that most models can only handle one response variable (default: %(default)s)')
    parser.add_argument('-ml', '--model-library', default="weka",
                        help='Model visitor library to use. (default: %(default)s)')
    parser.add_argument('-mt', '--model-tool', default="SVM_PUK_basic",
                        help='Model visitor tool from library to use. (default: %(default)s)')
    parser.add_argument('-s', '--splitter', default="KFoldSplitter",
                        help='Splitter to use. (default: %(default)s)')
    parser.add_argument('-v', dest='verbose', action='store_true',
                        help='Activate verbose mode.')
    args = parser.parse_args()

    prepare_build_display_model(args.predictor_headers, args.response_headers, args.model_library, args.model_tool, args.splitter, verbose=args.verbose)
