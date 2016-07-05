"""A module containing the ModelContainer class and related classes for descriptor attributes."""
from django.db import models
from django.conf import settings
from django.db import transaction
from django.core.exceptions import ValidationError
from numpy import average
from itertools import chain, zip_longest
import random
import datetime
import importlib
import os
from DRP.models.rxnDescriptors import BoolRxnDescriptor, OrdRxnDescriptor, NumRxnDescriptor, CatRxnDescriptor
from DRP.models.rxnDescriptorValues import BoolRxnDescriptorValue, NumRxnDescriptorValue, OrdRxnDescriptorValue, CatRxnDescriptorValue
from .statsModel import StatsModel
from DRP.utils import accuracy, BCR, Matthews, confusionMatrixString, confusionMatrixTable
import json
import sys
import logging

logger = logging.getLogger(__name__)

visitorModules = {library: importlib.import_module(
    settings.STATS_MODEL_LIBS_DIR + "." + library) for library in settings.STATS_MODEL_LIBS}

splitters = {splitter: importlib.import_module(
    settings.REACTION_DATASET_SPLITTERS_DIR + "." + splitter) for splitter in settings.REACTION_DATASET_SPLITTERS}

featureVisitorModules = {library: importlib.import_module(
    settings.FEATURE_SELECTION_LIBS_DIR + "." + library) for library in settings.FEATURE_SELECTION_LIBS}

MODEL_VISITOR_TOOL_CHOICES = tuple(
    tool for library in visitorModules.values() for tool in library.tools)

FEATURE_SELECTION_TOOL_CHOICES = tuple(
    tool for library in featureVisitorModules.values() for tool in library.tools)


class PredictsDescriptorsAttribute(object):
    """An attribute manager object which allows the setting and deletion of the related predictable descriptors."""

    def __get__(self, modelContainer, modelContainerType=None):
        """Return a chain of each of the querysets."""
        return chain(modelContainer.predboolrxndescriptor_set.all(), modelContainer.predordrxndescriptor_set.all(), modelContainer.predcatrxndescriptor_set.all(), modelContainer.prednumrxndescriptor_set.all())

    def __set__(self, modelContainer, descriptors):
        """Clear the descriptor sets and then set them equal to the provided queryset."""
        modelContainer.predboolrxndescriptor_set.clear()
        modelContainer.predordrxndescriptor_set.clear()
        modelContainer.predcatrxndescriptor_set.clear()
        modelContainer.prednumrxndescriptor_set.clear()
        for descriptor in descriptors:
            descriptor.stats_modelContainer = modelContainer
            descriptors.save()

    def __delete__(self, modelContainer):
        """Simply clear the querysets."""
        modelContainer.predboolrxndescriptor_set.clear()
        modelContainer.predordrxndescriptor_set.clear()
        modelContainer.predcatrxndescriptor_set.clear()
        modelContainer.prednumrxndescriptor_set.clear()


class DescriptorAttribute(object):
    """An attribute manager object which allows the setting and deletion of the related descriptors."""

    def __get__(self, modelContainer, modelContainerType=None):
        """Return a chain of each of the querysets."""
        return chain(modelContainer.boolRxnDescriptors.all(), modelContainer.ordRxnDescriptors.all(), modelContainer.catRxnDescriptors.all(), modelContainer.numRxnDescriptors.all())

    def __set__(self, modelContainer, descriptors):
        """Clear the descriptor sets and then set them equal to the provided queryset."""
        modelContainer.boolRxnDescriptors.clear()
        modelContainer.ordRxnDescriptors.clear()
        modelContainer.catRxnDescriptors.clear()
        modelContainer.numRxnDescriptors.clear()
        for descriptor in descriptors:
            # downcasting. Yuck.
            desc = None
            try:
                desc = BoolRxnDescriptor.objects.get(id=descriptor.id)
                modelContainer.boolRxnDescriptors.add(desc)
            except BoolRxnDescriptor.DoesNotExist:
                pass

            try:
                desc = OrdRxnDescriptor.objects.get(id=descriptor.id)
                modelContainer.ordRxnDescriptors.add(desc)
            except OrdRxnDescriptor.DoesNotExist:
                pass

            try:
                desc = CatRxnDescriptor.objects.get(id=descriptor.id)
                modelContainer.catRxnDescriptors.add(desc)
            except CatRxnDescriptor.DoesNotExist:
                pass

            try:
                desc = NumRxnDescriptor.objects.get(id=descriptor.id)
                modelContainer.numRxnDescriptors.add(desc)
            except NumRxnDescriptor.DoesNotExist:
                pass

            if desc is None:
                raise ValueError(
                    'An invalid object was assigned as a descriptor')

    def __delete__(self, modelContainer):
        """Simply clear the querysets."""
        modelContainer.boolRxnDescriptors.clear()
        modelContainer.numRxnDescriptors.clear()
        modelContainer.catRxnDescriptors.clear()
        modelContainer.ordRxnDescriptors.clear()


class OutcomeDescriptorAttribute(object):
    """An attribute manager object which allows the setting and deletion of the related outcome descriptors."""

    def __get__(self, modelContainer, modelContainerType=None):
        """Return a chain of each of the querysets."""
        return chain(modelContainer.outcomeBoolRxnDescriptors.all(), modelContainer.outcomeOrdRxnDescriptors.all(), modelContainer.outcomeCatRxnDescriptors.all(), modelContainer.outcomeNumRxnDescriptors.all())

    def __set__(self, modelContainer, descriptors):
        """Clear the descriptor sets and then set them equal to the provided queryset."""
        modelContainer.outcomeBoolRxnDescriptors.clear()
        modelContainer.outcomeOrdRxnDescriptors.clear()
        modelContainer.outcomeCatRxnDescriptors.clear()
        modelContainer.outcomeNumRxnDescriptors.clear()
        for descriptor in descriptors:
            desc = None
            try:
                desc = BoolRxnDescriptor.objects.get(id=descriptor.id)
                modelContainer.outcomeBoolRxnDescriptors.add(desc)
            except BoolRxnDescriptor.DoesNotExist:
                pass

            try:
                desc = OrdRxnDescriptor.objects.get(id=descriptor.id)
                modelContainer.outcomeOrdRxnDescriptors.add(desc)
            except OrdRxnDescriptor.DoesNotExist:
                pass

            try:
                desc = CatRxnDescriptor.objects.get(id=descriptor.id)
                modelContainer.outcomeCatRxnDescriptors.add(desc)
            except CatRxnDescriptor.DoesNotExist:
                pass

            try:
                desc = NumRxnDescriptor.objects.get(id=descriptor.id)
                modelContainer.outcomeNumRxnDescriptors.add(desc)
            except NumRxnDescriptor.DoesNotExist:
                pass

            if desc is None:
                raise ValueError(
                    'An invalid object was assigned as a descriptor')

    def __delete__(self, modelContainer):
        """Simply clear the querysets."""
        modelContainer.outcomeBoolRxnDescriptors.clear()
        modelContainer.outcomeNumRxnDescriptors.clear()
        modelContainer.outcomeCatRxnDescriptors.clear()
        modelContainer.outcomeOrdRxnDescriptors.clear()


class ModelContainer(models.Model):
    """A class for describing a group of statistical models."""

    class Meta:
        app_label = 'DRP'

    description = models.TextField(blank=True, null=False)
    active = models.BooleanField('Is this the active model?', default=False)
    modelVisitorLibrary = models.CharField(max_length=200)
    # choices=tuple((lib, lib) for lib in settings.STATS_MODEL_LIBS)
    modelVisitorTool = models.CharField(max_length=200)
    # choices=tuple((tool, tool) for tool in MODEL_VISITOR_TOOL_CHOICES)
    splitter = models.CharField(max_length=200, blank=True, default='')
    # choices=tuple((splitter, splitter) for splitter in settings.REACTION_DATASET_SPLITTERS)

    modelVisitorOptions = models.TextField(
        null=False, blank=True, default="{}")
    splitterOptions = models.TextField(null=False, blank=True, default="{}")

    built = models.BooleanField(
        'Has the build procedure been called with this container?', editable=False, default=False)

    descriptors = DescriptorAttribute()
    boolRxnDescriptors = models.ManyToManyField(BoolRxnDescriptor)
    ordRxnDescriptors = models.ManyToManyField(OrdRxnDescriptor)
    catRxnDescriptors = models.ManyToManyField(CatRxnDescriptor)
    numRxnDescriptors = models.ManyToManyField(NumRxnDescriptor)
    """The input descriptors for the model."""

    outcomeDescriptors = OutcomeDescriptorAttribute()
    outcomeBoolRxnDescriptors = models.ManyToManyField(
        BoolRxnDescriptor, related_name='outcomeForModels')
    outcomeOrdRxnDescriptors = models.ManyToManyField(
        OrdRxnDescriptor, related_name='outcomeForModels')
    outcomeCatRxnDescriptors = models.ManyToManyField(
        CatRxnDescriptor, related_name='outcomeForModels')
    outcomeNumRxnDescriptors = models.ManyToManyField(
        NumRxnDescriptor, related_name='outcomeForModels')
    """The descriptors which are being used as outcomes for this model.

    For models which make predictions about descriptors, it is probably
    most appropriate to make the descriptor label for the predicted
    descriptor related to the model, for instance for an outcome
    descriptor called "outcome", you might consider
    "outcome_predicted_by_model_id_1", where 1 is the model primary
    key.
    """

    predictsDescriptors = PredictsDescriptorsAttribute()
    """The descriptors which this model predicts values for."""

    fully_trained = models.ForeignKey("DRP.StatsModel", null=True, blank=True)

    @classmethod
    def create(cls, modelVisitorLibrary, modelVisitorTool, predictors, responses, splitterOptions=None, visitorOptions=None, description="",
               splitter=None, reactions=None, trainingSets=None, testSets=[], verbose=False):
        """
        Create a model container with a set of options.

        modelVisitorLibrary sepcifies the library/plugin being used, while the modelVisitorTool should specify the actual tool.
        predictors is the descriptors being used as independent variables, responses is a descriptor queryset for the independent variables.
        splitterOptions and visitorOptions depend on the splitter and vistior being used.
        the splitter is an imported class which instructs data how to be divided into training and test sets.
        reactions or training sets should be specified, and this will define how the training data is defined, as is testSets.

        """
        model_container = cls(modelVisitorLibrary=modelVisitorLibrary,
                              modelVisitorTool=modelVisitorTool, description=description)

        if (splitter is None) ^ (reactions is None):  # if these are not the same, there's a problem
            raise ValidationError(
                'A full set of reactions must be supplied with a splitter', 'argument_mismatch')
        # if these are not different, there's a problem
        if not ((splitter is None) ^ (trainingSets is None)):
            raise ValidationError(
                'Either a splitter or a training set should be provided.', 'argument_mismatch')

        if splitterOptions is not None:
            if splitter is None:
                raise ValidationError(
                    'Cannot define splitter options with no splitter')
            else:
                model_container.splitterOptions = json.dumps(splitterOptions)
        if visitorOptions is not None:
            model_container.modelVisitorOptions = json.dumps(visitorOptions)

        model_container.save()

        if splitter is not None:
            model_container.splitter = splitter
            splitter_name_stub = "{}_{}_{}".format(
                model_container.modelVisitorLibrary, model_container.modelVisitorTool, model_container.pk)
            splitterObj = splitters[model_container.splitter].Splitter(
                splitter_name_stub, **json.loads(model_container.splitterOptions))
            if verbose:
                logger.info("Splitting using {}".format(
                    model_container.splitter))
            data_splits = splitterObj.split(reactions, verbose=verbose)
        else:
            if verbose:
                logger.info("Using given test and training sets.")
            # we want the trainingset even if there's no test set
            data_splits = zip_longest(trainingSets, testSets)

        data_splits = model_container.createStatsModels(
            data_splits, verbose=verbose)

        model_container.descriptors = predictors
        model_container.outcomeDescriptors = responses

        return model_container

    def create_duplicate(self, modelVisitorTool=None, modelVisitorOptions=None, description=None, predictors=None, responses=None):
        """
        Build a duplicate of this model container optionally with a different model visitor tool.

        If a new description is not specified then the old description is used with 'rebuilt with tool X' appended
        """
        fields = ['description', 'splitter', 'splitterOptions',
                  'modelVisitorLibrary', 'modelVisitorTool', 'modelVisitorOptions']
        field_dict = ModelContainer.objects.filter(
            pk=self.pk).values(*fields)[0]

        if modelVisitorTool is not None:
            field_dict['modelVisitorTool'] = modelVisitorTool
        if modelVisitorOptions is not None:
            field_dict['modelVisitorOptions'] = json.dumps(modelVisitorOptions)
        if description is not None:
            field_dict['description'] = description
        else:
            addendum = u' rebuilt'
            if modelVisitorTool is not None:
                addendum += ' with tool {}'.format(modelVisitorTool)
            if modelVisitorOptions is not None:
                addendum += ' with options {}'.format(modelVisitorOptions)
            field_dict['description'] += addendum

        m = ModelContainer(**field_dict)
        m.save()

        if predictors is None:
            m.descriptors = self.descriptors
        else:
            m.descriptors = predictors
        if responses is None:
            m.outcomeDescriptors = self.outcomeDescriptors
        else:
            m.outcomeDescriptors = responses

        if self.statsmodel_set.all():
            for sm in self.statsmodel_set.all():
                statsModel = StatsModel(
                    container=m, trainingSet=sm.trainingSet)
                if set(m.descriptors) == set(self.descriptors) and set(m.outcomeDescriptors) == set(self.outcomeDescriptors):
                    statsModel.inputFile = sm.inputFile
                statsModel.save()
                statsModel.testSets = sm.testSets.all()
        else:
            raise RuntimeError(
                'This model container was never properly constructed, so it cannot be duplicated. (It has no stats models)')

        return m

    def clean(self):
        """Perform very rudimentary validation of additional properties. Needs refactoring."""
        if self.modelVisitorTool not in visitorModules[self.modelVisitorLibrary].tools:
            raise ValidationError('Selected tool {} does not exist in selected library {}'.format(
                self.modelVisitorTool, self.modelVisitorLibrary), 'wrong_library')
        if getattr(visitorModules[self.modelVisitorLibrary], self.modelVisitorTool).maxResponseCount is not None:
            if getattr(visitorModules[self.modelVisitorLibrary], self.modelVisitorTool).maxResponseCount < len([d for d in self.outcomeDescriptors]):
                raise ValidationError('Selected tool {} cannot accept this many responses, maximum is {}'.format(self.modelVisitorTool, getattr(
                    visitorModules[self.modelVisitorLibrary], self.modelVisitorTool).maxResponseCount), 'too_many_responses')
        try:
            options_dict = json.loads(self.modelVisitorOptions)
        except:
            raise ValidationError('Was unable to parse modelVisitorOptions {} with json. Got exception: ({})'.format(
                self.modelVisitorOptions, repr(sys.exc_info()[1])))
        try:
            getattr(visitorModules[self.modelVisitorLibrary], self.modelVisitorTool)(
                statsModel=None, **options_dict)
        except:
            raise ValidationError('Was unable expand modelVisitorOptions {} parsed by json into keyword arguments accepted by model visitor. Got exception: {}'.format(
                self.modelVisitorOptions, repr(sys.exc_info()[1])))
        try:
            options_dict = json.loads(self.splitterOptions)
        except:
            raise ValidationError('Was unable to parse splitterOptions {} with json. Got exception: ({})'.format(
                self.splitterOptions, repr(sys.exc_info()[1])))
        try:
            splitterObj = splitters[self.splitter].Splitter('', **options_dict)
        except:
            raise ValidationError('Was unable to expand splitterOptions {} parsed by json into keyword arguments accepted by splitter. Got exception: {}'.format(
                self.splitterOptions, repr(sys.exc_info()[1])))

    def createStatsModels(self, data_splits, verbose=False):
        """Create statistical models which 'vote' on the outcomes. May be singular or multiple."""
        for trainingSet, testSet in data_splits:
            statsModel = StatsModel(container=self, trainingSet=trainingSet)
            statsModel.save()
            statsModel.testSets.add(testSet)

    def build(self, verbose=False):
        """
        Take all options confirmed so far and generate a full model set.

        Train a mutlitude of models using the external libraries selected, then save the
        relevant information to the database (see the statsmodel class.)

        Run the tests for the model using the test sets of data, and then saves that information.

        """
        if self.built:
            raise RuntimeError(
                "Cannot build a model that has already been built.")

        if verbose:
            logger.info("Starting building at {}".format(
                datetime.datetime.now()))

        # set up a prediction results dictionary. Hold on tight. This gets
        # hairy real fast.
        resDict = {}

        num_models = self.statsmodel_set.all().count()
        num_finished = 0
        overall_start_time = datetime.datetime.now()
        for statsModel in self.statsmodel_set.all():
            visitorOptions = json.loads(self.modelVisitorOptions)
            modelVisitor = getattr(visitorModules[self.modelVisitorLibrary], self.modelVisitorTool)(
                statsModel=statsModel, **visitorOptions)
            # Train the model.
            statsModel.startTime = datetime.datetime.now()
            fileName = os.path.join(settings.MODEL_DIR, '{}_{}_{}_{}.model'.format(
                self.pk, statsModel.pk, self.modelVisitorLibrary, self.modelVisitorTool))
            statsModel.outputFile = fileName
            if verbose:
                logger.info("{} statsModel {}, saving to {}, training...".format(
                    statsModel.startTime, statsModel.pk, fileName))
            modelVisitor.train(verbose=verbose)
            statsModel.endTime = datetime.datetime.now()
            if verbose:
                logger.info("\t...Trained. Finished at {}. Saving statsModel...".format(
                    statsModel.endTime)),
            statsModel.save()
            if verbose:
                logger.info("saved")

            # Test the model.
            for testSet in statsModel.testSets.all():
                if testSet.reactions.all().count() != 0:
                    if verbose:
                        logger.info("Predicting test set...")
                    predictions = modelVisitor.predict(
                        testSet.reactions.all(), verbose=verbose)
                    if verbose:
                        logger.info(
                            "\t...finished predicting. Storing predictions...",)
                    newResDict = self._storePredictionComponents(
                        predictions, statsModel)

                    # Update the overall result-dictionary with these new
                    # counts.
                    for reaction, responseDict in newResDict.items():
                        for response, outcomeDict in responseDict.items():
                            for outcome, count in outcomeDict.items():
                                if reaction not in resDict:
                                    resDict[reaction] = {}
                                if response not in resDict[reaction]:
                                    resDict[reaction][response] = {}

                                if outcome not in resDict[reaction][response]:
                                    resDict[reaction][response][
                                        outcome] = count
                                resDict[reaction][response][outcome] += count

                    if verbose:
                        logger.info("predictions stored.")
                        for response in self.outcomeDescriptors:
                            predDesc = response.predictedDescriptorType.objects.get(
                                modelContainer=self, statsModel=statsModel, predictionOf=response)
                            conf_mtrx = predDesc.getConfusionMatrix()

                            logger.info(
                                "Confusion matrix for {}:".format(predDesc.heading))
                            logger.info(confusionMatrixString(conf_mtrx))
                            logger.info("Accuracy: {:.3}".format(
                                accuracy(conf_mtrx)))
                            logger.info("BCR: {:.3}".format(BCR(conf_mtrx)))

                elif verbose:
                    logger.info("Test set is empty.")

            if verbose:
                num_finished += 1
                end_time = datetime.datetime.now()
                elapsed = (end_time - overall_start_time)
                expected_finish = datetime.timedelta(seconds=(elapsed.total_seconds(
                ) * (num_models / float(num_finished)))) + overall_start_time
                logger.info("{}. {} of {} models built.".format(
                    end_time, num_finished, num_models))
                logger.info("Elapsed model building time: {}. Expected completion time: {}".format(
                    elapsed, expected_finish))

        if resDict:
            if verbose:
                logger.info("Storing overall model predictions...")
            self._storePredictions(resDict)
            if verbose:
                logger.info("Predictions stored")

        self.built = True
        if verbose:
            overall_end_time = datetime.datetime.now()
            logger.info("Finished at {}".format(overall_end_time))

    def _storePredictionComponents(self, predictions, statsModel, resDict=None):
        """
        Return resDict, a dictionary of dictionaries of dictionaries.

        The first key is the reaction, the second is the response descriptor
        (the descriptor to be predicted), the third is the predicted outcome.
        The value is 1 if that outcome is predicted and 0 otherwise.
        This setup is for easier aggregating into the full model container
        voting-based prediction.
        """
        resDict = {} if resDict is None else resDict

        for response, outcomes in predictions.items():
            predDesc = response.createPredictionDescriptor(self, statsModel)
            predDesc.save()
            vals = []
            for reaction, outcome in outcomes:
                if reaction not in resDict:
                    resDict[reaction] = {}
                if response not in resDict[reaction]:
                    resDict[reaction][response] = {}
                if outcome not in resDict[reaction][response]:
                    resDict[reaction][response][outcome] = 0
                resDict[reaction][response][outcome] += 1
                val = predDesc.createValue(reaction, outcome)
                if val.pk is None:  # if the value already exists
                    vals.append(val)
                else:
                    val.save()
            if vals:
                type(vals[0]).objects.bulk_create(vals)
        return resDict

    def _storePredictions(self, resDict):
        """Store predictions from the overall container as voted for by each componenet model."""
        finalPredictions = {}
        bool_vals = []
        num_vals = []
        ord_vals = []
        cat_vals = []
        for reaction, responseDict in resDict.items():
            for response, outcomeDict in responseDict.items():
                predDesc = response.createPredictionDescriptor(self)
                if predDesc.pk is None:
                    predDesc.save()
                if isinstance(response, NumRxnDescriptor):
                    # estimate the weighted average of the estimates from the
                    # component models
                    values = tuple(value for value in resDict[
                                   reaction][response].keys())
                    weights = tuple(weight for value, weight in resDict[
                                    reaction][response].items())
                    val = predDesc.createValue(
                        reaction, average(values, weights=weights))
                    num_vals.append(val)
                else:
                    # find the 'competitors' with the highest number of votes
                    # and pick one at random (in case multiple categories have
                    # equal numbers of votes)
                    maxVotes = max(count for response,
                                   count in outcomeDict.items())
                    winners = [(response, count) for response,
                               count in outcomeDict.items() if count == maxVotes]
                    winner = random.choice(winners)
                    val = predDesc.createValue(reaction, winner[0])
                    if val.pk is not None:
                        val.save()
                    elif isinstance(val, BoolRxnDescriptorValue):
                        bool_vals.append(val)
                    elif isinstance(val, OrdRxnDescriptorValue):
                        ord_vals.append(val)
                    elif isinstance(val, CatRxnDescriptorValue):
                        cat_vals.append(val)
                    else:
                        raise ValueError(
                            "Value just created is of unexpected type {}".format(type(val)))

                if response not in finalPredictions:
                    finalPredictions[response] = []
                finalPredictions[response].append((reaction, val))

        BoolRxnDescriptorValue.objects.bulk_create(bool_vals)
        NumRxnDescriptorValue.objects.bulk_create(num_vals)
        OrdRxnDescriptorValue.objects.bulk_create(ord_vals)
        CatRxnDescriptorValue.objects.bulk_create(cat_vals)
        return finalPredictions

    def predict(self, reactions, verbose=False):
        """Make predictions from the voting for a set of provided reactions."""
        if self.built:
            resDict = {}

            num_models = self.statsmodel_set.all().count()
            num_finished = 0
            overall_start_time = datetime.datetime.now()
            for model in self.statsmodel_set.all():
                visitorOptions = json.loads(self.modelVisitorOptions)
                modelVisitor = getattr(visitorModules[self.modelVisitorLibrary], self.modelVisitorTool)(
                    statsModel=model, **visitorOptions)
                if verbose:
                    logger.info("statsModel {}, saved at {}, predicting...".format(
                        model.pk, model.outputFile))
                predictions = modelVisitor.predict(reactions, verbose=verbose)
                if verbose:
                    logger.info(
                        "\t...finished predicting. Storing predictions...")
                newResDict = self._storePredictionComponents(
                    predictions, model)

                # Update the overall result-dictionary with these new counts.
                for reaction, responseDict in newResDict.items():
                    for response, outcomeDict in responseDict.items():
                        for outcome, count in outcomeDict.items():
                            if reaction not in resDict:
                                resDict[reaction] = {}
                            if response not in resDict[reaction]:
                                resDict[reaction][response] = {}

                            if outcome not in resDict[reaction][response]:
                                resDict[reaction][response][outcome] = count
                            resDict[reaction][response][outcome] += count

                if verbose:
                    logger.info("predictions stored.")
                    for response in self.outcomeDescriptors:
                        predDesc = response.predictedDescriptorType.objects.get(
                            modelContainer=self, statsModel=model, predictionOf=response)
                        conf_mtrx = predDesc.getConfusionMatrix(
                            reactions=reactions)

                        logger.info(
                            "Confusion matrix for {}:".format(predDesc.heading))
                        logger.info(confusionMatrixString(conf_mtrx))
                        logger.info("Accuracy: {:.3}".format(
                            accuracy(conf_mtrx)))
                        logger.info("BCR: {:.3}".format(BCR(conf_mtrx)))
                        logger.info("Matthews: {:.3}".format(
                            Matthews(conf_mtrx)))
                    num_finished += 1
                    end_time = datetime.datetime.now()
                    elapsed = (end_time - overall_start_time)
                    expected_finish = datetime.timedelta(seconds=(elapsed.total_seconds(
                    ) * (num_models / float(num_finished)))) + overall_start_time
                    logger.info("{}. Predictions from {} of {} models.".format(
                        end_time, num_finished, num_models))
                    logger.info("Elapsed prediction time: {}. Expected completion time: {}".format(
                        elapsed, expected_finish))

            return self._storePredictions(resDict)
        else:
            raise RuntimeError(
                'A model container cannot be used to make predictions before the build method has been called')

    def getOverallConfusionMatrices(self, reactions=None):
        """Return the confusion matrix for the voted predictions from this ModelContainer."""
        confusion_matrix_list = []
        for descriptor in self.predictsDescriptors:
            if descriptor.statsModel is None:
                confusion_matrix_list.append(
                    (descriptor.csvHeader, descriptor.getConfusionMatrix(reactions=reactions)))
        return confusion_matrix_list

    def getComponentConfusionMatrices(self, reactions=None):
        """
        Return a list of lists of tuples of confusion matrices.

        Each entry of the outer list is for a different component statsModel.
        For each model there is a list of tuples.
        Each tuple is of the form (descriptor_heading, confusion matrix)
        """
        confusion_matrix_lol = []

        for model in self.statsmodel_set.all():
            confusion_matrix_list = []
            for descriptor in self.predictsDescriptors:
                if descriptor.statsModel == model:
                    confusion_matrix_list.append(
                        (descriptor.csvHeader, descriptor.getConfusionMatrix(reactions=reactions)))
            confusion_matrix_lol.append(confusion_matrix_list)

        return confusion_matrix_lol
