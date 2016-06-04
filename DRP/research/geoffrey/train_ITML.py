#!/usr/bin/env python

import django
import argparse
from DRP.research.geoffrey.distance_learning.metricLearn import ITML
from DRP.models import PerformedReaction, ModelContainer, Descriptor, rxnDescriptorValues

def train(predictor_headers, response_headers, outfile, num_constraints):
    # Grab all valid reactions with defined outcome descriptors
    reactions = PerformedReaction.objects.filter(valid=True)
    reactions = reactions.exclude(ordrxndescriptorvalue__in=rxnDescriptorValues.OrdRxnDescriptorValue.objects.filter(descriptor__heading__in=response_headers, value=None))
    reactions = reactions.exclude(boolrxndescriptorvalue__in=rxnDescriptorValues.BoolRxnDescriptorValue.objects.filter(descriptor__heading__in=response_headers, value=None))
    reactions = reactions.exclude(catrxndescriptorvalue__in=rxnDescriptorValues.CatRxnDescriptorValue.objects.filter(descriptor__heading__in=response_headers, value=None))

    predictors = Descriptor.objects.filter(heading__in=predictor_headers)
    responses = Descriptor.objects.filter(heading__in=response_headers)

    predictor_headers = [d.csvHeader for d in predictors]
    response_headers = [d.csvHeader for d in responses]
    
    itml = ITML()

    itml.train(reactions, predictor_headers, response_headers, num_constraints=num_constraints)

    with open(outfile, 'wb') as f:
        itml.save(f)


if __name__ == '__main__':
    django.setup()
    parser = argparse.ArgumentParser(description='Trains an ITML metric', fromfile_prefix_chars='@',
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
    parser.add_argument('-d', '--description', default="",
                        help='Description of model. (default: %(default)s)')
    parser.add_argument('-o', '--outfile', required=True,
                        help='Place to dump the pickle.')
    parser.add_argument('-n', '--num-constraints', required=True, type=int,
                        help='Number of constraints for the ITML')
    args = parser.parse_args()

    train(args.predictor_headers, args.response_headers, args.outfile, args.num_constraints)
