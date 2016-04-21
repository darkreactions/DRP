#!/usr/bin/env python

import django
from DRP.models import PerformedReaction, ModelContainer, Descriptor, rxnDescriptorValues, DataSet
import operator
import argparse
from django.db.utils import OperationalError
#from DRP.ml_models.splitters.MutualInfoSplitter import Splitter
import uuid
from itertools import chain
from django.conf import settings
import importlib
import ast

splitters = {splitter:importlib.import_module(settings.REACTION_DATASET_SPLITTERS_DIR + "." + splitter) for splitter in settings.REACTION_DATASET_SPLITTERS}

def prepareArff(reactions, whitelistHeaders, description, verbose=False):
    """Writes an *.arff file using the provided queryset of reactions."""
    filepath = "{}_{}.arff".format(description, uuid.uuid4())
    if verbose:
        print "Writing arff to {}".format(filepath)
    with open(filepath, "w") as f:
        reactions.toArff(f, expanded=True, whitelistHeaders=whitelistHeaders)
    return filepath

def split_and_dump(predictor_headers=None, response_headers=None, reaction_set_name=None, description="", verbose=False, splitterOptions={}, splitter=None):
    predictors = Descriptor.objects.filter(heading__in=predictor_headers)
    responses = Descriptor.objects.filter(heading__in=response_headers)

    if predictors.count() != len(predictor_headers):
        raise KeyError("Could not find all predictors")
    if responses.count() != len(response_headers):
        raise KeyError("Could not find all responses")

    if reaction_set_name is not None:
        reactions = DataSet.objects.get(name=reaction_set_name).reactions.all()
    else:
        reactions = PerformedReaction.objects.all()

    splitter_name_stub = "{}_{}".format(description, uuid.uuid4())
    splitterObj = splitters[splitter].Splitter(splitter_name_stub, **splitterOptions)
    if verbose:
        print "Splitting using {}".format(splitter)
    data_splits = splitterObj.split(reactions, verbose=verbose)

    whitelist = [d.csvHeader for d in chain(predictors, responses)]

    for trainingSet, testSet in data_splits:
        prepareArff(trainingSet.reactions.all(), whitelist, trainingSet.name + '_train', verbose=verbose)
        prepareArff(testSet.reactions.all(), whitelist, testSet.name + '_test', verbose=verbose)

if __name__ == '__main__':
    django.setup()
    parser = argparse.ArgumentParser(description='Splits reaction set into datasets and dumps to arffs', fromfile_prefix_chars='@',
                                     epilog="Prefix arguments with '@' to specify a file containing newline"
                                     "-separated values for that argument. e.g.'-p @predictor_headers.txt'"
                                     " to pass multiple descriptors from a file as predictors")
    parser.add_argument('-p', '--predictor-headers', nargs='+',
                        help='One or more descriptors to use as predictors.', required=True)
    parser.add_argument('-r', '--response-headers', nargs='+', default=["boolean_crystallisation_outcome"],
                        help='One or more descriptors to predict. '
                        'Note that most models can only handle one response variable (default: %(default)s)')
    parser.add_argument('-s', '--splitter', default="KFoldSplitter",
                        help='Splitter to use. (default: %(default)s)')
    parser.add_argument('-so', '--splitter-options', default=None,
                        help='A dictionary of the options to give to the splitter in JSON format')
    parser.add_argument('-v', dest='verbose', action='store_true',
                        help='Activate verbose mode.')
    parser.add_argument('-d', '--description', default="",
                        help='Description of dataset. (default: %(default)s)')
    parser.add_argument('-rxn', '--reaction-set-name', default=None,
                        help='The name of the reactions to use as a whole dataset')
    args = parser.parse_args()

    splitterOptions = ast.literal_eval(args.splitter_options) if args.splitter_options is not None else {}

    split_and_dump(predictor_headers=args.predictor_headers, response_headers=args.response_headers, reaction_set_name=args.reaction_set_name,
                description=args.description, verbose=args.verbose, splitterOptions=splitterOptions, splitter=args.splitter)
