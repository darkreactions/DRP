"""Command for building statistical/machine learning models in DRP."""
from django.core.management.base import BaseCommand
from DRP.models import PerformedReaction, ModelContainer, Descriptor, rxnDescriptorValues, DataSet
import operator
import argparse
from django.db.utils import OperationalError
from time import sleep
from DRP.utils import accuracy, BCR, Matthews, confusionMatrixString, confusionMatrixTable
from django.conf import settings
import ast
from sys import argv
import logging
logger = logging.getLogger(__name__)


class Command(BaseCommand):
    """Command for building statistical/machine learning models in DRP."""

    help = 'Builds a machine learning model'

    def add_arguments(self, parser):
        """Add arguments for the argument parser."""
        parser.fromfile_prefix_chars = '@'
        parser.epilog = ("Prefix arguments with '@' to specify a file containing newline-separated values for that argument. "
                         "e.g.'-p @predictor_headers.txt' to pass multiple descriptor headings from a file as predictors.")

        parser.add_argument('-p', '--predictor-headers', nargs='+',
                            help='The headings of one or more descriptors to use as predictors.')
        parser.add_argument('-r', '--response-headers', nargs='+', default=["boolean_crystallisation_outcome"],
                            help='The headings of one or more descriptors to predict. '
                            'Note that most models can only handle one response variable (default: %(default)s)')
        parser.add_argument('-ml', '--model-library', default="weka",
                            help='Model visitor library to use. (default: %(default)s)')
        parser.add_argument('-mt', '--model-tool', default="SVM_PUK",
                            help='Model visitor tool from library to use. (default: %(default)s)')
        parser.add_argument('-mid', '--model-container-id', default=None, type=int,
                            help='Use the same splits as the specified model container. (default: %(default)s)')
        parser.add_argument('-s', '--splitter', default="KFoldSplitter", choices=settings.REACTION_DATASET_SPLITTERS,
                            help='Splitter to use. (default: %(default)s)')
        parser.add_argument('-d', '--description', default="",
                            help='Description of model. (default: %(default)s)')
        parser.add_argument('-trs', '--training-set-name', default=None,
                            help='The name of the training set to use. (default: %(default)s)')
        parser.add_argument('-tes', '--test-set-name', default=None,
                            help='The name of the test set to use. (default: %(default)s)')
        parser.add_argument('-rxn', '--reaction-set-name', default=None,
                            help='The name of the reactions to use as a whole dataset. (default: all valid reactions)')
        parser.add_argument('-so', '--splitter-options', default=None,
                            help='A dictionary of the options to give to the splitter in JSON format')
        parser.add_argument('-vo', '--visitor-options', default=None,
                            help='A dictionary of the options to give to the visitor in JSON format')

        # TODO setup argparse to properly check combinations of these arguments as valid.
        # it's actually pretty complicated what the valid options are...

    def handle(self, *args, **kwargs):
        """Handle the call for this command."""
        # This way of accepting splitter options is bad and hacky.
        # Unfortunately, the only good ways I can think of are also very complicated and I don't have time right now :-(
        # TODO XXX make this not horrible
        splitterOptions = ast.literal_eval(kwargs['splitter_options']) if kwargs[
            'splitter_options'] is not None else None
        visitorOptions = ast.literal_eval(kwargs['visitor_options']) if kwargs[
            'visitor_options'] is not None else None

        # Remove errant empty strings
        predictor_headers = [h for h in kwargs['predictor_headers'] if h] if kwargs[
            'predictor_headers'] is not None else None
        response_headers = [h for h in kwargs['response_headers'] if h] if kwargs[
            'response_headers'] is not None else None

        # TODO switch to logging and adjust for management command multi-level
        # verbosity
        verbose = (kwargs['verbosity'] > 0)

        prepare_build_display_model(predictor_headers=predictor_headers, response_headers=response_headers,
                                    modelVisitorLibrary=kwargs[
                                        'model_library'], modelVisitorTool=kwargs['model_tool'],
                                    splitter=kwargs['splitter'], training_set_name=kwargs['training_set_name'], test_set_name=kwargs['test_set_name'], reaction_set_name=kwargs['reaction_set_name'], description=kwargs['description'], verbose=verbose, container_id=kwargs['model_container_id'], splitterOptions=splitterOptions, visitorOptions=visitorOptions)


# TODO refactor this so these functions can be used in other places more naturally
# Also so this workflow makes some freaking sense. What the hell was I thinking?
# Handle should call each part in turn, not have a bunch of parameters that keep getting passed through till someone uses them
# I think I did that so that other files could call model building at the point it was need and jump into the pipeline,
# but wow there has to be a better way - GCMN

def create_build_model(reactions=None, predictors=None, responses=None, modelVisitorLibrary=None, modelVisitorTool=None, splitter=None, trainingSet=None, testSet=None,
                       description=None, verbose=False, splitterOptions=None, visitorOptions=None):
    """Build the model and puts it into the DB."""
    if trainingSet is not None:
        container = ModelContainer.create(modelVisitorLibrary, modelVisitorTool, predictors, responses, description=description, reactions=reactions,
                                          trainingSets=[trainingSet], testSets=[
                                              testSet], verbose=verbose, splitterOptions=splitterOptions,
                                          visitorOptions=visitorOptions)
    else:
        container = ModelContainer.create(modelVisitorLibrary, modelVisitorTool, predictors, responses, description=description, reactions=reactions,
                                          splitter=splitter, verbose=verbose, splitterOptions=splitterOptions, visitorOptions=visitorOptions)

    container.full_clean()
    return build_model(container, verbose=verbose)


def build_model(container, verbose=False):
    """An additional function by GMN to build models. I don't really know what it's for- PA."""
    for attempt in range(5):
        try:
            container.build(verbose=verbose)
            break
        except OperationalError as e:
            logger.warning("Caught OperationalError {}\nRestarting in 3 seconds...\n".format(e))
            sleep(3)
    else:
        raise RuntimeError("Got 5 Operational Errors in a row and gave up")

    container.save()
    container.full_clean()

    return container


def missing_descriptors(descriptor_headings):
    """Find descriptors which were requested but aren't in the DB.'."""
    missing_descs = []
    for heading in descriptor_headings:
        if not Descriptor.objects.filter(heading=heading).exists():
            missing_descs.append(heading)
    return missing_descs


def display_model_results(container, reactions=None, heading=""):
    """
    Display confusion matrices for a model container.

    Optional heading specifies prefix for the summary statistics
    (useful for when multiple model containers are built by a single script)
    """
    overall_conf_mtrcs = container.getOverallConfusionMatrices(
        reactions=reactions)
    if not overall_conf_mtrcs:
        print("No model results to display")
        return
    if len(overall_conf_mtrcs) != 1:
        raise NotImplementedError('Can only handle one response')
    for descriptor_header, conf_mtrx in overall_conf_mtrcs:
        acc = accuracy(conf_mtrx)
        bcr = BCR(conf_mtrx)
        matthews = Matthews(conf_mtrx)
        print("Confusion matrix for {}:".format(descriptor_header))
        print(confusionMatrixString(conf_mtrx))
        print("Accuracy: {:.3}".format(acc))
        print("BCR: {:.3}".format(bcr))
        print("Matthews: {:.3}".format(matthews))
    conf_mtrcs = container.getComponentConfusionMatrices(reactions=reactions)

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
            print("Confusion matrix for {}:".format(descriptor_header))
            print(confusionMatrixString(conf_mtrx))
            print("Accuracy: {:.3}".format(acc))
            print("BCR: {:.3}".format(bcr))
            print("Matthews: {:.3}".format(matthews))

            # This only works for one response. Sorry...
            # TODO XXX make this work for multiple responses
            sum_acc += acc
            sum_bcr += bcr
            sum_matthews += matthews
            count += 1

    print("{} Average accuracy: {:.3}".format(heading, sum_acc / count))
    print("{} Average BCR: {:.3}".format(heading, sum_bcr / count))
    print("{} Average Matthews: {:.3}".format(heading, sum_matthews / count))


def prepare_build_model(predictor_headers=None, response_headers=None, modelVisitorLibrary=None, modelVisitorTool=None, splitter=None, training_set_name=None,
                        test_set_name=None, reaction_set_name=None, description=None, verbose=False, splitterOptions=None, visitorOptions=None, container_id=None):
    """Build a model with the specified tools."""
    if predictor_headers is not None:
        predictors = Descriptor.objects.filter(heading__in=predictor_headers)
        if predictors.count() != len(predictor_headers):
            raise KeyError("Could not find all predictors. Missing: {}".format(
                missing_descriptors(predictor_headers)))
    else:
        predictors = None
    if response_headers is not None:
        responses = Descriptor.objects.filter(heading__in=response_headers)
        if responses.count() != len(response_headers):
            raise KeyError("Could not find all responses. Missing: {}".format(
                missing_descriptors(response_headers)))
    else:
        responses = None

    if container_id is not None:
        parent_container = ModelContainer.objects.get(id=container_id)
        parent_container.full_clean()
        new_container = parent_container.create_duplicate(
            modelVisitorTool=modelVisitorTool, modelVisitorOptions=visitorOptions, description=description, predictors=predictors, responses=responses)
        new_container.full_clean()
        container = build_model(new_container, verbose=verbose)
    else:
        if training_set_name is None and reaction_set_name is None:
            assert(test_set_name is None)
            reactions = PerformedReaction.objects.filter(valid=True)
            trainingSet = None
            testSet = None
        elif reaction_set_name is not None:
            reaction_set = DataSet.objects.get(name=reaction_set_name)
            reactions = reaction_set.reactions.all()
            trainingSet = None
            testSet = None
        else:
            trainingSet = DataSet.objects.get(name=training_set_name)
            testSet = DataSet.objects.get(name=test_set_name)
            reactions = None

        container = create_build_model(reactions=reactions, predictors=predictors, responses=responses,
                                       modelVisitorLibrary=modelVisitorLibrary, modelVisitorTool=modelVisitorTool,
                                       splitter=splitter, trainingSet=trainingSet, testSet=testSet,
                                       description=description, verbose=verbose, splitterOptions=splitterOptions,
                                       visitorOptions=visitorOptions)

    return container


def prepare_build_display_model(predictor_headers=None, response_headers=None, modelVisitorLibrary=None, modelVisitorTool=None, splitter=None, training_set_name=None, test_set_name=None,
                                reaction_set_name=None, description=None, verbose=False, splitterOptions=None, visitorOptions=None, container_id=None):
    """I'm not exactly clear on what this function by GMN is for- PA."""
    container = prepare_build_model(predictor_headers=predictor_headers, response_headers=response_headers, modelVisitorLibrary=modelVisitorLibrary, modelVisitorTool=modelVisitorTool,
                                    splitter=splitter, training_set_name=training_set_name, test_set_name=test_set_name, reaction_set_name=reaction_set_name, description=description,
                                    verbose=verbose, splitterOptions=splitterOptions, visitorOptions=visitorOptions, container_id=container_id)

    display_model_results(container)
