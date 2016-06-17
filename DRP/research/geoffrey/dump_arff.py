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
import build_model


def prepareArff(reactions, whitelistHeaders, filepath, verbose=False):
    """Writes an *.arff file using the provided queryset of reactions."""
    if verbose:
        print "Writing arff to {}".format(filepath)
    with open(filepath, "w") as f:
        reactions.toArff(f, expanded=True, whitelistHeaders=whitelistHeaders)
    return filepath


def dump(response_headers=None, predictor_headers=None, reaction_set_name=None, description="", verbose=False, output_file=None):
    responses = Descriptor.objects.filter(heading__in=response_headers)
    predictors = Descriptor.objects.filter(heading__in=predictor_headers)

    if responses.count() != len(response_headers):
        missing_descriptors = build_model.missing_descriptors(response_headers)
        raise KeyError(
            "Could not find all responses. Missing {}".format(missing_descriptors))
    if predictors.count() != len(predictor_headers):
        missing_descriptors = build_model.missing_descriptors(
            predictor_headers)
        raise KeyError(
            "Could not find all predictors. Missing {}".format(missing_descriptors))

    if reaction_set_name is None:
        reactions = PerformedReaction.objects.all()
    else:
        reactions = DataSet.objects.get(name=reaction_set_name).reactions.all()

    whitelist = [d.csvHeader for d in chain(predictors, responses)]
    prepareArff(reactions, whitelist, output_file, verbose=verbose)


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
    parser.add_argument('-v', dest='verbose', action='store_true',
                        help='Activate verbose mode.')
    parser.add_argument('-rxn', '--reaction-set-name', default=None,
                        help='The name of the reactions to use as a whole dataset')
    parser.add_argument('-o', '--output-file', default=None, required=True,
                        help='The file to dump the arff to')

    args = parser.parse_args()

    dump(predictor_headers=args.predictor_headers, response_headers=args.response_headers,
         reaction_set_name=args.reaction_set_name, output_file=args.output_file, verbose=args.verbose)
