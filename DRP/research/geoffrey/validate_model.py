#!/usr/bin/env python

import django
from DRP.models import PerformedReaction, ModelContainer, Descriptor, DataSet, NumRxnDescriptor, BoolRxnDescriptor
import argparse
from django.conf import settings
import ast
from sys import argv
import json
import build_model
from django.db.models import Count

def check_container(container, predictor_headers=None, response_headers=None, modelVisitorLibrary=None, modelVisitorTool=None, splitter=None,
                                verbose=False, splitterOptions=None, visitorOptions=None):
    raise NotImplementedError("Haven't implemented checking a given container")

def find_container(model_id=None, predictor_headers=None, response_headers=None, modelVisitorLibrary=None, modelVisitorTool=None, splitter=None,
                                verbose=False, unique=False, splitterOptions=None, visitorOptions=None):

    if model_id is not None:
        try:
            container = ModelContainer.objects.get(id=model_id)
            check_container(container, predictor_headers, response_headers, modelVisitorLibrary, modelVisitorTool, splitter, verbose, splitterOptions, visitorOptions)
        except ModelContainer.DoesNotExist:
            raise RuntimeError("No model container with that id exists")
    else:
        containers = ModelContainer.objects.filter(built=True)
        if verbose:
            print "Found {} built model containers".format(containers.count())

        if splitter is not None:
            containers = containers.filter(splitter=splitter)
            if verbose:
                print "Found {} with correct splitter".format(containers.count())
        if splitterOptions is not None:
            containers = containers.filter(splitterOptions=json.dumps(splitterOptions))
            if verbose:
                print "Found {} with correct splitter options".format(containers.count())


        if modelVisitorLibrary is not None:
            containers = containers.filter(modelVisitorLibrary=modelVisitorLibrary)
            if verbose:
                print "Found {} with correct visitor library".format(containers.count())
        if modelVisitorTool is not None:
            containers = containers.filter(modelVisitorTool=modelVisitorTool)
            if verbose:
                print "Found {} with correct visitor tool".format(containers.count())
        if visitorOptions is not None:
            containers = containers.filter(modelVisitorOptions=json.dumps(visitorOptions))
            if verbose:
                print "Found {} with correct visitor options".format(containers.count())

        if response_headers is not None:
            containers = containers.annotate(num_outNumDescs=Count('outcomeNumRxnDescriptors', distinct=True)).annotate(num_outBoolDescs=Count('outcomeBoolRxnDescriptors', distinct=True))
            boolDescs = BoolRxnDescriptor.objects.filter(heading__in=response_headers)
            numDescs = NumRxnDescriptor.objects.filter(heading__in=response_headers)
            num_numDescs = numDescs.count()
            num_boolDescs = boolDescs.count()
    
            setOutBoolDescs = set(boolDescs)
            setOutNumDescs = set(numDescs)

            if num_numDescs + num_boolDescs != len(response_headers):
                missing_descs = []
                for header in headers:
                    if not NumRxnDescriptor.objects.filter(heading=header).exists() and not BoolRxnDescriptor.objects.filter(heading=header).exists():
                        missing_descs.append(header)
                raise RuntimeError("Did not find correct number of descriptors. Currently only looking at boolean and numeric. Categorical and ordinal unimplemented. Unable to find:{}".format(missing_descs))
    
            containers = containers.filter(num_outNumDescs=num_numDescs).filter(num_outBoolDescs=num_boolDescs)
            if verbose:
                print "{} containers with appropriate number of responses".format(containers.count())

        if predictor_headers is not None:
            containers = containers.annotate(num_numDescs=Count('numRxnDescriptors', distinct=True)).annotate(num_boolDescs=Count('boolRxnDescriptors', distinct=True))
            boolDescs = BoolRxnDescriptor.objects.filter(heading__in=predictor_headers)
            numDescs = NumRxnDescriptor.objects.filter(heading__in=predictor_headers)
            num_numDescs = numDescs.count()
            num_boolDescs = boolDescs.count()
    
            setBoolDescs = set(boolDescs)
            setNumDescs = set(numDescs)

            if num_numDescs + num_boolDescs != len(predictor_headers):
                missing_descs = []
                for header in headers:
                    if not NumRxnDescriptor.objects.filter(heading=header).exists() and not BoolRxnDescriptor.objects.filter(heading=header).exists():
                        missing_descs.append(header)
                raise RuntimeError("Did not find correct number of descriptors. Currently only looking at boolean and numeric. Categorical and ordinal unimplemented. Unable to find:{}".format(missing_descs))
    
            containers = containers.filter(num_numDescs=num_numDescs).filter(num_boolDescs=num_boolDescs)
            if verbose:
                print "{} containers with appropriate number of predictors".format(containers.count())

        if response_headers is not None:
            containers = [c for c in containers if (set(c.outcomeNumRxnDescriptors.all()) == setOutNumDescs and set(c.outcomeBoolRxnDescriptors.all()) == setOutBoolDescs)]
            if verbose:
                print "{} containers with correct responses".format(len(containers))
            
            
        if predictor_headers is not None:
            containers = [c for c in containers if (set(c.numRxnDescriptors.all()) == setNumDescs and set(c.boolRxnDescriptors.all()) == setBoolDescs)]
            if verbose:
                print "{} containers with correct predictors".format(len(containers))

        if len(containers) != 1:
            print "Was unable to find a unique model container matching given specification. Found {}".format(len(containers))
            if len(containers) > 1 and not unique:
                print "Using container with largest pk"
                try:
                    containers.sort(key=lambda x: x.pk, reverse=True)
                except AttributeError:
                    containers = containers.order_by('-pk')
            else:
                raise RuntimeError("No unique container. If you want a non-unique container, remove the '-u' flag")

        return containers[0]

def find_predict_display_model(model_id=None, predictor_headers=None, response_headers=None, modelVisitorLibrary=None, modelVisitorTool=None, splitter=None, reaction_set_name=None,
                                verbose=False, unique=False, splitterOptions=None, visitorOptions=None):

    container = find_container(model_id=model_id, predictor_headers=predictor_headers, response_headers=response_headers, modelVisitorLibrary=modelVisitorLibrary,
                                modelVisitorTool=modelVisitorTool, splitter=splitter, verbose=verbose, unique=unique, splitterOptions=splitterOptions, visitorOptions=visitorOptions)

    reaction_set = DataSet.objects.get(name=reaction_set_name)
    reactions = reaction_set.reactions.all()

    container.predict(reactions)

    build_model.display_model_results(container,reactions)


if __name__ == '__main__':
    django.setup()
    parser = argparse.ArgumentParser(description='Test a model on a give reaction set', fromfile_prefix_chars='@',
                                     epilog="Prefix arguments with '@' to specify a file containing newline"
                                     "-separated values for that argument. e.g.'-p @predictor_headers.txt'"
                                     " to pass multiple descriptors from a file as predictors")
    parser.add_argument('-id', '--model-id')
    parser.add_argument('-p', '--predictor-headers', nargs='+',
                        help='One or more descriptors to use as predictors.')
    parser.add_argument('-r', '--response-headers', nargs='+', default=["boolean_crystallisation_outcome"],
                        help='One or more descriptors to predict. '
                        'Note that most models can only handle one response variable (default: %(default)s)')
    parser.add_argument('-ml', '--model-library', default="weka",
                        help='Model visitor library to use. (default: %(default)s)')
    parser.add_argument('-mt', '--model-tool', default="SVM_PUK",
                        help='Model visitor tool from library to use. (default: %(default)s)')
    parser.add_argument('-s', '--splitter', default="KFoldSplitter", choices=settings.REACTION_DATASET_SPLITTERS,
                        help='Splitter to use. (default: %(default)s)')
    parser.add_argument('-v', '--verbose', action='store_true',
                        help='Activate verbose mode.')
    parser.add_argument('-u', '--unique', action='store_true',
                        help='Find only unique container.')
    parser.add_argument('-so', '--splitter-options', default=None,
                        help='A dictionary of the options to give to the splitter in JSON format')
    parser.add_argument('-vo', '--visitor-options', default=None,
                        help='A dictionary of the options to give to the visitor in JSON format')
    parser.add_argument('-rxn', '--reaction-set-name', default=None, required=True,
                        help='The name of the reaction set to predict')
                        
    args = parser.parse_args()
    if args.verbose:
        print argv[1:]
        print args


    # This way of accepting splitter options is bad and hacky.
    # Unfortunately, the only good ways I can think of are also very complicated and I don't have time right now :-(
    # TODO XXX make this not horrible
    splitterOptions = ast.literal_eval(args.splitter_options) if args.splitter_options is not None else None
    visitorOptions = ast.literal_eval(args.visitor_options) if args.visitor_options is not None else None

    find_predict_display_model(model_id=args.model_id, predictor_headers=args.predictor_headers, response_headers=args.response_headers, modelVisitorLibrary=args.model_library, modelVisitorTool=args.model_tool,
                                splitter=args.splitter, verbose=args.verbose, splitterOptions=splitterOptions, visitorOptions=visitorOptions, unique=args.unique, reaction_set_name=args.reaction_set_name)
