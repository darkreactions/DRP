#!/usr/bin/env python

import django
from DRP.models import PerformedReaction, MetricContainer, Descriptor, rxnDescriptorValues, DataSet
import argparse
from DRP.ml_models.splitters.SingleSplitter import Splitter as SingleSplitter

def prepare_build_metric(descriptor_headers=None, response_headers=None, metricVisitorTool=None, description="", trainingSetName=None, testSetName=None, outfile=None, verbose=False):
    """
    Build and display a model with the specified tools
    """

    trainingSet = DataSet.objects.get(name=trainingSetName)
    predictors = Descriptor.objects.filter(heading__in=descriptor_headers)
    responses = Descriptor.objects.filter(heading__in=response_headers)

    container = MetricContainer(metricVisitor=metricVisitorTool, trainingSet=trainingSet, description=description)
    container.save()
    container.full_clean()
    transformed = container.build(predictors, responses, verbose=verbose)
    container.save()
    container.full_clean()


    if verbose:
        print "Transforming training set"
    reactions = trainingSet.reactions.all()
    container.transform(reactions, transformed=transformed, verbose=verbose)

    if verbose:
        print "Writing descriptors to file"

    if outfile is not None:
        with open(outfile, 'wb') as f:
            for desc in container.transformedRxnDescriptors.all():
                f.write(desc.heading)
                f.write('\n')

    if testSetName is not None:
        if verbose:
            print "Tranforming test set to new space" 
        testSet = DataSet.objects.get(name=testSetName)
        reactions = testSet.reactions.all()
        container.transform(reactions, verbose=verbose)


if __name__ == '__main__':
    django.setup()
    parser = argparse.ArgumentParser(description='Builds a model', fromfile_prefix_chars='@',
                                     epilog="Prefix arguments with '@' to specify a file containing newline"
                                     "-separated values for that argument. e.g.'-p @predictor_headers.txt'"
                                     " to pass multiple descriptors from a file as predictors")
    parser.add_argument('-p', '--predictor-headers', nargs='+',
                        help='One or more descriptors to use as predictors.', required=True)
    parser.add_argument('-r', '--response-headers', nargs='+', default=["boolean_crystallisation_outcome"],
                        help='One or more descriptors to predict. '
                        'Note that most models can only handle one response variable (default: %(default)s)')
    parser.add_argument('-m', '--metric-tool', default="ITML",
                        help='Metric visitor to use. (default: %(default)s)')
    parser.add_argument('-v', dest='verbose', action='store_true',
                        help='Activate verbose mode.')
    parser.add_argument('-d', '--description', default="",
                        help='Description of metric.')
    parser.add_argument('-o', '--descriptor-outfile', default=None,
                        help='File to write list of metric descriptors.')
    parser.add_argument('-trs', '--training-set-name', default=None, required=True,
                        help='Name of training set to use.')
    parser.add_argument('-tes', '--test-set-name', default=None, required=True,
                        help='Name of training set to use.')
    args = parser.parse_args()

    prepare_build_metric(args.predictor_headers, args.response_headers, args.metric_tool, args.description, trainingSetName=args.training_set_name, testSetName=args.test_set_name, outfile=args.descriptor_outfile, verbose=args.verbose)
