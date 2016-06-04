#!/usr/bin/env python

import django
from DRP.models import PerformedReaction, ModelContainer, Descriptor, rxnDescriptorValues, DataSet
import operator
import argparse
from django.db.utils import OperationalError
from DRP.ml_models.splitters.MutualInfoSplitter import Splitter
import uuid
from itertools import chain
import glob
import os


def prepareArff(reactions, whitelistHeaders, filepath, verbose=False):
    """Writes an *.arff file using the provided queryset of reactions."""
    if verbose:
        print "Writing arff to {}".format(filepath)
    with open(filepath, "w") as f:
        reactions.toArff(f, expanded=True, whitelistHeaders=whitelistHeaders)
    return filepath


def split_and_dump(response_headers=None, reaction_set_name=None, description="", verbose=False):
    responses = Descriptor.objects.filter(heading__in=response_headers)

    if responses.count() != len(response_headers):
        raise KeyError("Could not find all responses")

    if reaction_set_name is not None:
        reactions = DataSet.objects.get(name=reaction_set_name).reactions.all()
    else:
        reactions = PerformedReaction.objects.all()
        reaction_set_name = "all"

    splitterObj = Splitter("{}_{}".format(reaction_set_name, uuid.uuid4()))

    data_splits = splitterObj.split(reactions, verbose=verbose)

    filepaths = sorted(glob.glob('../legacy_tests/descs/step*.dsc'))

    with open('../legacy_tests/descs/new_full.dsc') as f:
        old_predictor_headers = [l.strip() for l in f.readlines() if l.strip()]

    old_predictors = Descriptor.objects.filter(heading__in=old_predictor_headers)

    if old_predictors.count() != len(old_predictor_headers):
        raise KeyError("Could not find all predictors")

    for trainingSet, testSet in data_splits:
        for desc_file in filepaths:
            with open(desc_file) as f:
                predictor_headers = [l.strip() for l in f.readlines() if l.strip()]
            predictors = Descriptor.objects.filter(heading__in=predictor_headers)

            if predictors.count() != len(predictor_headers):
                raise KeyError("Could not find all predictors")

            whitelist = [d.csvHeader for d in chain(old_predictors, predictors, responses)]
            desc_file_description = os.path.basename(desc_file)[:-4]
            prepareArff(trainingSet.reactions.all(), whitelist, "{}_{}_train.arff".format(desc_file_description, trainingSet.name), verbose=verbose)
            prepareArff(testSet.reactions.all(), whitelist, "{}_{}_test.arff".format(desc_file_description, testSet.name), verbose=verbose)

if __name__ == '__main__':
    django.setup()
    parser = argparse.ArgumentParser(description='Splits reaction set into datasets and dumps to arffs', fromfile_prefix_chars='@',
                                     epilog="Prefix arguments with '@' to specify a file containing newline"
                                     "-separated values for that argument. e.g.'-p @predictor_headers.txt'"
                                     " to pass multiple descriptors from a file as predictors")
    parser.add_argument('-r', '--response-headers', nargs='+', default=["boolean_crystallisation_outcome"],
                        help='One or more descriptors to predict. '
                        'Note that most models can only handle one response variable (default: %(default)s)')
    parser.add_argument('-v', dest='verbose', action='store_true',
                        help='Activate verbose mode.')
    parser.add_argument('-rxn', '--reaction-set-name', default=None,
                        help='The name of the reactions to use as a whole dataset')
    args = parser.parse_args()

    split_and_dump(response_headers=args.response_headers, reaction_set_name=args.reaction_set_name, verbose=args.verbose)
