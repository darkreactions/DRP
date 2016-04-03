#!/usr/bin/env python

import django
from DRP.models import PerformedReaction, ModelContainer, Descriptor, rxnDescriptorValues, DataSet
import argparse
import build_model
from django.conf import settings
from itertools import izip
import ast


def prepare_build_display_many_models(predictor_headers=None, response_headers=None, modelVisitorLibrary=None, modelVisitorTools=None, splitter=None, training_set_name=None, test_set_name=None,
                                      reaction_set_name=None, description="", verbose=False, splitter_options=None, visitor_options=None):

    if visitor_options is not None and len(modelVisitorTools) != len(visitor_options):
        raise ValidationError('Need to specify options for all models or none')

    visitors_with_options = izip(modelVisitorTools, visitor_options)
    
    initialVisitor, initial_visitor_options = visitors_with_options.next()
    
    container = build_model.prepare_build_model(predictor_headers=predictor_headers, response_headers=response_headers, modelVisitorLibrary=modelVisitorLibrary, modelVisitorTool=initialVisitor,
                                    splitter=splitter, training_set_name=training_set_name, test_set_name=test_set_name, reaction_set_name=reaction_set_name, description=description,
                                    verbose=verbose, splitter_options=splitter_options, visitor_options=initial_visitor_options)

    build_model.display_model_results(container)

    for visitor, options in visitors_with_options:
        visitorOptions = ast.literal_eval(options)

        new_container = container.create_duplicate(modelVisitorTool=visitor, modelVisitorOptions=visitorOptions)
        new_container.full_clean()

        build_model.build_model(new_container, verbose=verbose)
        build_model.display_model_results(new_container)
        
    

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
    parser.add_argument('-ml', '--model-library', default="weka",
                        help='Model visitor library to use. (default: %(default)s)')
    parser.add_argument('-mt', '--model-tools', nargs='+',
                        help='Model visitor tools from library to use.')
    parser.add_argument('-s', '--splitter', default="KFoldSplitter", choices=settings.REACTION_DATASET_SPLITTERS,
                        help='Splitter to use. (default: %(default)s)')
    parser.add_argument('-v', '--verbose', action='store_true',
                        help='Activate verbose mode.')
    parser.add_argument('-d', '--description', default="",
                        help='Description of model. (default: %(default)s)')
    parser.add_argument('-trs', '--training-set-name', default=None,
                        help='The name of the training set to use. (default: %(default)s)')
    parser.add_argument('-tes', '--test-set-name', default=None,
                        help='The name of the test set to use. (default: %(default)s)')
    parser.add_argument('-rxn', '--reaction-set-name', default=None,
                        help='The name of the reactions to use as a whole dataset')
    parser.add_argument('-so', '--splitter-options', default=None,
                        help='A dictionary of the options to give to the splitter in JSON format')
    parser.add_argument('-vo', '--visitor-options', default=None, nargs='+',
                        help='A dictionary of the options to give to each visitor in JSON format')
                        
    args = parser.parse_args()

    prepare_build_display_many_models(predictor_headers=args.predictor_headers, response_headers=args.response_headers, modelVisitorLibrary=args.model_library, modelVisitorTools=args.model_tools,
                                splitter=args.splitter, training_set_name=args.training_set_name, test_set_name=args.test_set_name, reaction_set_name=args.reaction_set_name,
                                description=args.description, verbose=args.verbose, splitter_options=args.splitter_options, visitor_options=args.visitor_options)
