#!/usr/bin/env python

import django
from DRP.models import PerformedReaction, ModelContainer, Descriptor, rxnDescriptorValues, DataSet
import operator
import argparse
from django.db.utils import OperationalError
from time import sleep
from DRP.utils import accuracy, BCR, Matthews, confusionMatrixString, confusionMatrixTable
from django.conf import settings
import ast
from sys import argv

def create_build_model(reactions=None, predictors=None, responses=None, modelVisitorLibrary=None, modelVisitorTool=None, splitter=None, trainingSet=None, testSet=None,
                description="", verbose=False, splitterOptions=None, visitorOptions=None):

    if trainingSet is not None:
        container = ModelContainer.create(modelVisitorLibrary, modelVisitorTool, predictors, responses, description=description, reactions=reactions,
                                          trainingSets=[trainingSet], testSets=[testSet], verbose=verbose, splitterOptions=splitterOptions,
                                          visitorOptions=visitorOptions)
    else:
        container = ModelContainer.create(modelVisitorLibrary, modelVisitorTool, predictors, responses, description=description, reactions=reactions,
                                          splitter=splitter, verbose=verbose, splitterOptions=splitterOptions, visitorOptions=visitorOptions)

    container.full_clean()
    return build_model(container, verbose=verbose)

def build_model(container, verbose=False):
    for attempt in range(5):
        try:
            container.build(verbose=verbose)
            break
        except OperationalError as e:
            print "Caught OperationalError {}".format(e)
            print "\nRestarting in 3 seconds...\n"
            sleep(3)
    else:
        raise RuntimeError("Got 5 Operational Errors in a row and gave up")
    
    container.save()
    container.full_clean()

    return container

def missing_descriptors(descriptor_headings):
    missing_descs = []
    for heading in descriptor_headings:
        if not Descriptor.objects.filter(heading=heading).exists():
            missing_descs.append(heading)
    return missing_descs

def display_model_results(container, heading=""):
    """
    Displays confusion matrices for a model container.
    Optional heading specifies prefix for the summary statistics
    (useful for when multiple model containers are built by a single script)
    """
    overall_conf_mtrcs = container.getOverallConfusionMatrices()
    if not overall_conf_mtrcs:
        print "No model results to display"
        return
    if len(overall_conf_mtrcs) != 1:
        raise NotImplementedError('Can only handle one response')
    for descriptor_header, conf_mtrx in overall_conf_mtrcs:
        acc = accuracy(conf_mtrx)
        bcr = BCR(conf_mtrx)
        matthews = Matthews(conf_mtrx)
        print "Confusion matrix for {}:".format(descriptor_header)
        print confusionMatrixString(conf_mtrx)
        print "Accuracy: {:.3}".format(acc)
        print "BCR: {:.3}".format(bcr)
        print "Matthews: {:.3}".format(matthews)
    conf_mtrcs = container.getComponentConfusionMatrices()

    sum_acc = 0.0
    sum_bcr = 0.0
    sum_matthews = 0.0
    count = 0

    for model_mtrcs in conf_mtrcs:
        if len(model_mtrcs) != 1:
            raise NotImplementedError('Can only handle one response')
        for descriptor_header, conf_mtrx in model_mtrcs:
            acc = accuracy(conf_mtrx)
            bcr = BCR(conf_mtrx)
            matthews = Matthews(conf_mtrx)
            print "Confusion matrix for {}:".format(descriptor_header)
            print confusionMatrixString(conf_mtrx)
            print "Accuracy: {:.3}".format(acc)
            print "BCR: {:.3}".format(bcr)
            print "Matthews: {:.3}".format(matthews)

            # This only works for one response. Sorry...
            # TODO XXX make this work for multiple responses
            sum_acc += acc
            sum_bcr += bcr
            sum_matthews += matthews
            count += 1
            
    print "{} Average accuracy: {:.3}".format(heading, sum_acc/count)
    print "{} Average BCR: {:.3}".format(heading, sum_bcr/count)
    print "{} Average Matthews: {:.3}".format(heading, sum_matthews/count)

def prepare_build_model(predictor_headers=None, response_headers=None, modelVisitorLibrary=None, modelVisitorTool=None, splitter=None, training_set_name=None,
                        test_set_name=None, reaction_set_name=None, description="", verbose=False, splitterOptions=None, visitorOptions=None):
    """
    Build a model with the specified tools
    """

    # Remove errant empty strings
    predictor_headers = [h for h in predictor_headers if h]
    response_headers = [h for h in response_headers if h]
    
    predictors = Descriptor.objects.filter(heading__in=predictor_headers)
    responses = Descriptor.objects.filter(heading__in=response_headers)

    if predictors.count() != len(predictor_headers):
        raise KeyError("Could not find all predictors. Missing: {}".format(missing_descriptors(predictor_headers)))
    if responses.count() != len(response_headers):
        raise KeyError("Could not find all responses. Missing: {}".format(missing_descriptors(response_headers)))

    if training_set_name is None and reaction_set_name is None:
        assert(test_set_name == None)
        reactions = PerformedReaction.objects.filter(valid=True)
        #reactions = reactions.exclude(ordrxndescriptorvalue__in=rxnDescriptorValues.OrdRxnDescriptorValue.objects.filter(descriptor__heading__in=response_headers, value=None))
        #reactions = reactions.exclude(boolrxndescriptorvalue__in=rxnDescriptorValues.BoolRxnDescriptorValue.objects.filter(descriptor__heading__in=response_headers, value=None))
        #reactions = reactions.exclude(catrxndescriptorvalue__in=rxnDescriptorValues.CatRxnDescriptorValue.objects.filter(descriptor__heading__in=response_headers, value=None))
        trainingSet = None
        testSet = None
    elif reaction_set_name is not None:
        reaction_set = DataSet.objects.get(name=reaction_set_name)
        reactions = reaction_set.reactions.all()
        trainingSet = None
        testSet = None
    else:
        trainingSet =  DataSet.objects.get(name=training_set_name)
        testSet =  DataSet.objects.get(name=test_set_name)
        reactions=None
    
    container = create_build_model(reactions=reactions, predictors=predictors, responses=responses, 
                            modelVisitorLibrary=modelVisitorLibrary, modelVisitorTool=modelVisitorTool,
                            splitter=splitter, trainingSet=trainingSet, testSet=testSet, 
                            description=description, verbose=verbose, splitterOptions=splitterOptions,
                            visitorOptions=visitorOptions)

    return container
    
def prepare_build_display_model(predictor_headers=None, response_headers=None, modelVisitorLibrary=None, modelVisitorTool=None, splitter=None, training_set_name=None, test_set_name=None,
                                reaction_set_name=None, description="", verbose=False, splitterOptions=None, visitorOptions=None):

    container = prepare_build_model(predictor_headers=predictor_headers, response_headers=response_headers, modelVisitorLibrary=modelVisitorLibrary, modelVisitorTool=modelVisitorTool,
                                    splitter=splitter, training_set_name=training_set_name, test_set_name=test_set_name, reaction_set_name=reaction_set_name, description=description,
                                    verbose=verbose, splitterOptions=splitterOptions, visitorOptions=visitorOptions)

    display_model_results(container)


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
    parser.add_argument('-mt', '--model-tool', default="SVM_PUK",
                        help='Model visitor tool from library to use. (default: %(default)s)')
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
                        help='A dictionary of the options to give to the visitor in JSON format')
                        
    args = parser.parse_args()
    if args.verbose:
        print argv[1:]
        print args


    # This way of accepting splitter options is bad and hacky.
    # Unfortunately, the only good ways I can think of are also very complicated and I don't have time right now :-(
    # TODO XXX make this not horrible
    splitterOptions = ast.literal_eval(args.splitter_options) if args.splitter_options is not None else None
    visitorOptions = ast.literal_eval(args.visitor_options) if args.visitor_options is not None else None

    prepare_build_display_model(predictor_headers=args.predictor_headers, response_headers=args.response_headers, modelVisitorLibrary=args.model_library, modelVisitorTool=args.model_tool,
                                splitter=args.splitter, training_set_name=args.training_set_name, test_set_name=args.test_set_name, reaction_set_name=args.reaction_set_name, description=args.description, verbose=args.verbose, splitterOptions=splitterOptions, visitorOptions=visitorOptions)
