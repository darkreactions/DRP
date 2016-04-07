#!/usr/bin/env python

import django
from DRP.models import PerformedReaction, ModelContainer, Descriptor, rxnDescriptorValues, DataSet
import argparse
import build_model
from django.conf import settings
from itertools import izip
import ast
import sys

def prepare_build_display_many_models(predictor_headers=None, response_headers=None, splitter=None, training_set_name=None, test_set_name=None, puk_sigma=None, puk_omega=None,
                                      reaction_set_name=None, description="", verbose=False, splitterOptions=None, visitorOptions=None):

    if visitorOptions is None:
        visitorOptions = {}

    if 'puk_sigma' in visitorOptions or 'puk_omega' in visitorOptions:
        raise ValidationError('Do not specify PUK sigma or omega in visitor options. Instead specify minimum, maximum, and step to scan over')

    modelVisitorLibrary = 'weka'
    modelVisitorTool = 'SVM_PUK'

    sigma_min, sigma_max, sigma_step = puk_sigma
    omega_min, omega_max, omega_step = puk_omega

    if sigma_max < sigma_min or omega_max < omega_min:
        raise ValidationError("Sigma max and omega max must be greater than respective min")

    sigma = sigma_min
    omega = omega_min
    
    visitorOptions['puk_sigma'] = sigma
    visitorOptions['puk_omega'] = omega
    
    initialDescription = "{} {} {}".format(modelVisitorTool, visitorOptions, description)

    if verbose:
        print "Building initial container with {} {}".format(modelVisitorTool, visitorOptions)
    container = build_model.prepare_build_model(predictor_headers=predictor_headers, response_headers=response_headers, modelVisitorLibrary=modelVisitorLibrary, modelVisitorTool=modelVisitorTool,
                                    splitter=splitter, training_set_name=training_set_name, test_set_name=test_set_name, reaction_set_name=reaction_set_name, description=initialDescription,
                                    verbose=verbose, splitterOptions=splitterOptions, visitorOptions=visitorOptions)

    build_model.display_model_results(container, heading='sigma={} omega={}'.format(sigma, omega))

    omega += omega_step

    while sigma <= sigma_max:
        while omega <= omega_max:
            if verbose:
                print "Building container with sigma={} omega={}".format(sigma, omega)

            visitorOptions['puk_sigma'] = sigma
            visitorOptions['puk_omega'] = omega

            new_description = "{} {} {}".format(modelVisitorTool, visitorOptions, description)

            new_container = container.create_duplicate(modelVisitorOptions=visitorOptions, description=new_description)
            new_container.full_clean()
    
            build_model.build_model(new_container, verbose=verbose)
            build_model.display_model_results(new_container, heading='sigma={} omega={}'.format(sigma, omega))

            omega += omega_step
        omega = omega_min
        sigma += sigma_step
        

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
    #parser.add_argument('-ml', '--model-library', default="weka",
                        #help='Model visitor library to use. (default: %(default)s)')
    #parser.add_argument('-mt', '--model-tools', nargs='+',
                        #help='Model visitor tools from library to use.')
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
    parser.add_argument('-vo', '--visitor-options', default=None,
                        help='A dictionary of the options to give to every visitor in JSON format.'
                             ' This should not include the PUK parameters.')
    parser.add_argument('-ps', '--puk-sigma', default=None, nargs=3, required=True, type=float,
                        help='Specify min, max, and step for the PUK sigma parameter')
    parser.add_argument('-po', '--puk-omega', default=None, nargs=3, required=True, type=float,
                        help='Specify min, max, and step for the PUK omega parameter')
    
    args = parser.parse_args()
    if args.verbose:
        print sys.argv[1:]
        print args

    # This way of accepting splitter options is bad and hacky.
    # Unfortunately, the only good ways I can think of are also very complicated and I don't have time right now :-(
    # TODO XXX make this not horrible
    splitterOptions = ast.literal_eval(args.splitter_options) if args.splitter_options is not None else None
    visitorOptions = ast.literal_eval(args.visitor_options) if args.visitor_options is not None else None

    prepare_build_display_many_models(predictor_headers=args.predictor_headers, response_headers=args.response_headers, puk_sigma=args.puk_sigma, puk_omega=args.puk_omega,
                                splitter=args.splitter, training_set_name=args.training_set_name, test_set_name=args.test_set_name, reaction_set_name=args.reaction_set_name,
                                description=args.description, verbose=args.verbose, splitterOptions=splitterOptions, visitorOptions=visitorOptions)
