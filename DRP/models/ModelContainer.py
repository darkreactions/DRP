"""A module containing only the ModelContainer class."""
from django.db import models
from django.conf import settings
from django.db import transaction
from django.core.exceptions import ValidationError
from numpy import average
from itertools import chain
import random
import datetime
import importlib
import os
from DRP.models.rxnDescriptors import BoolRxnDescriptor, OrdRxnDescriptor, NumRxnDescriptor, CatRxnDescriptor
from DRP.models.rxnDescriptorValues import BoolRxnDescriptorValue, NumRxnDescriptorValue, OrdRxnDescriptorValue, CatRxnDescriptorValue
from StatsModel import StatsModel
from DRP.research.geoffrey.utils import accuracy, BCR, confusionMatrixString, confusionMatrixTable

visitorModules = {library:importlib.import_module(settings.STATS_MODEL_LIBS_DIR + "."+ library) for library in settings.STATS_MODEL_LIBS}

splitters = {splitter:importlib.import_module(settings.REACTION_DATASET_SPLITTERS_DIR + "." + splitter) for splitter in settings.REACTION_DATASET_SPLITTERS}
#TODO: set availability of manual splitting up

featureVisitorModules = {library:importlib.import_module(settings.FEATURE_SELECTION_LIBS_DIR + "." + library) for library in settings.FEATURE_SELECTION_LIBS}

MODEL_VISITOR_TOOL_CHOICES = tuple(tool for library in visitorModules.values() for tool in library.tools)

FEATURE_SELECTION_TOOL_CHOICES = tuple(tool for library in featureVisitorModules.values() for tool in library.tools)

class PredictsDescriptorsAttribute(object):

    def __get__(self, modelContainer, modelContainerType=None):
        return chain(modelContainer.predboolrxndescriptor_set.all(), modelContainer.predordrxndescriptor_set.all(), modelContainer.predcatrxndescriptor_set.all(), modelContainer.prednumrxndescriptor_set.all())

    def __set__(self, modelContainer, descriptors):
        modelContainer.predboolrxndescriptor_set.clear()
        modelContainer.predordrxndescriptor_set.clear()
        modelContainer.predcatrxndescriptor_set.clear()
        modelContainer.prednumrxndescriptor_set.clear()
        for descriptor in descriptors:
            descriptor.stats_modelContainer = modelContainer
            descriptors.save()

    def __delete__(self, modelContainer):
        modelContainer.predboolrxndescriptor_set.clear()
        modelContainer.predordrxndescriptor_set.clear()
        modelContainer.predcatrxndescriptor_set.clear()
        modelContainer.prednumrxndescriptor_set.clear()

class DescriptorAttribute(object):

    def __get__(self, modelContainer, modelContainerType=None):
        return chain(modelContainer.boolRxnDescriptors.all(), modelContainer.ordRxnDescriptors.all(), modelContainer.catRxnDescriptors.all(), modelContainer.numRxnDescriptors.all())

    def __set__(self, modelContainer, descriptors):
        modelContainer.boolRxnDescriptors.clear()
        modelContainer.ordRxnDescriptors.clear()
        modelContainer.catRxnDescriptors.clear()
        modelContainer.numRxnDescriptors.clear()
        for descriptor in descriptors:
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
                raise ValueError('An invalid object was assigned as a descriptor')

    def __delete__(self, modelContainer):
        modelContainer.boolRxnDescriptors.clear()
        modelContainer.numRxnDescriptors.clear()
        modelContainer.catRxnDescriptors.clear()
        modelContainer.ordRxnDescriptors.clear()

class OutcomeDescriptorAttribute(object):

    def __get__(self, modelContainer, modelContainerType=None):
        return chain(modelContainer.outcomeBoolRxnDescriptors.all(), modelContainer.outcomeOrdRxnDescriptors.all(), modelContainer.outcomeCatRxnDescriptors.all(), modelContainer.outcomeNumRxnDescriptors.all())

    def __set__(self, modelContainer, descriptors):
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
                raise ValueError('An invalid object was assigned as a descriptor')

    def __delete__(self, modelContainer):
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
    #choices=tuple((splitter, splitter) for splitter in settings.REACTION_DATASET_SPLITTERS)
    
    # TODO XXX these should be validated as json or implemented another way (e.g. key-value store in another table)
    #modelVisitorOptions = models.TextField(null=False, blank=True, default="")
    #splitterOptions = models.TextField(null=False, blank=True, default="")
    
    built = models.BooleanField('Has the build procedure been called with this container?', editable=False, default=False)

    descriptors = DescriptorAttribute()
    boolRxnDescriptors = models.ManyToManyField(BoolRxnDescriptor)
    ordRxnDescriptors = models.ManyToManyField(OrdRxnDescriptor)
    catRxnDescriptors = models.ManyToManyField(CatRxnDescriptor)
    numRxnDescriptors = models.ManyToManyField(NumRxnDescriptor)
    """The input descriptors for the model."""

    outcomeDescriptors = OutcomeDescriptorAttribute()
    outcomeBoolRxnDescriptors = models.ManyToManyField(BoolRxnDescriptor, related_name='outcomeForModels')
    outcomeOrdRxnDescriptors = models.ManyToManyField(OrdRxnDescriptor, related_name='outcomeForModels')
    outcomeCatRxnDescriptors = models.ManyToManyField(CatRxnDescriptor, related_name='outcomeForModels')
    outcomeNumRxnDescriptors = models.ManyToManyField(NumRxnDescriptor, related_name='outcomeForModels')
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
    def create(cls, modelVisitorLibrary, modelVisitorTool, description="", splitter=None, reactions=None, trainingSets=None, testSets=None):
        model_container = cls(modelVisitorLibrary=modelVisitorLibrary, modelVisitorTool=modelVisitorTool, description=description)

        if (splitter is None) ^ (reactions is None): # if these are not the same, there's a problem
            raise ValidationError('A full set of reactions must be supplied with a splitter', 'argument_mismatch')
        if not ((splitter is None) ^ (trainingSets is None)): # if these are not different, there's a problem
            raise ValidationError('Either a splitter or a training set should be provided.', 'argument_mismatch')

        if splitter is not None:
            model_container.splitter = splitter

        model_container.reactions = reactions
        model_container.trainingSets = trainingSets
        model_container.testSets = testSets

        return model_container

    def clean(self):
        if self.modelVisitorTool not in visitorModules[self.modelVisitorLibrary].tools:
            raise ValidationError('Selected tool does not exist in selected library', 'wrong_library')
        if getattr(visitorModules[self.modelVisitorLibrary], self.modelVisitorTool).maxResponseCount is not None:
            if getattr(visitorModules[self.modelVisitorLibrary], self.modelVisitorTool).maxResponseCount < len([d for d in self.outcomeDescriptors]):
                raise ValidationError('Selected tool cannot accept this many responses, maximum is {}', 'too_many_responses', tuple(visitorModules[self.modelVisitorLibrary], self.modelVisitorTool).maxResponseCount)

    def build(self, predictors, response, verbose=False):
        if self.built:
            raise RuntimeError("Cannot build a model that has already been built.")

        if verbose:
            print "Starting building at {}".format(datetime.datetime.now())

        self.descriptors = predictors
        self.outcomeDescriptors = response

        if self.splitter:
            splitterObj = splitters[self.splitter].Splitter("{}_{}_{}".format(self.modelVisitorLibrary, self.modelVisitorTool, self.pk))
            if verbose:
                print "Splitting using {}".format(self.splitter)
            data_splits = splitterObj.split(self.reactions, verbose=verbose)

            self.trainingSets = [dataset_tuple[0] for dataset_tuple in data_splits]
            self.testSets = [dataset_tuple[1] for dataset_tuple in data_splits]
        elif verbose:
            print "No splitter. Using given test and training sets."
        
        resDict = {} #set up a prediction results dictionary. Hold on tight. This gets hairy real fast.


        num_models = len(self.trainingSets)
        num_finished = 0
        overall_start_time = datetime.datetime.now()
        for trainingSet, testSet in zip(self.trainingSets, self.testSets):
            statsModel = StatsModel(container=self, trainingSet=trainingSet)
            modelVisitor = getattr(visitorModules[self.modelVisitorLibrary], self.modelVisitorTool)(statsModel)
            statsModel.save() #Generate a PK for the StatsModel component.

            # Train the model.
            statsModel.startTime = datetime.datetime.now()
            fileName = os.path.join(settings.MODEL_DIR, '{}_{}_{}_{}.model'.format(self.pk, statsModel.pk, self.modelVisitorLibrary, self.modelVisitorTool))
            whitelist = [d.csvHeader for d in chain(self.descriptors, self.outcomeDescriptors)]
            if verbose:
                print "{} statsModel {}, saving to {}, training...".format(statsModel.startTime, statsModel.pk, fileName)
            modelVisitor.train(trainingSet.reactions.all(), whitelist, fileName, verbose=verbose)
            statsModel.outputFile = fileName
            statsModel.endTime = datetime.datetime.now()
            if verbose:
                print "\t...Trained. Finished at {}. Saving statsModel...".format(statsModel.endTime),
            statsModel.save()
            if verbose:
                print "saved"
            

            # Test the model.
            if testSet.reactions.count() > 0:
                if verbose:
                    print "Predicting test set..."
                statsModel.testSets.add(testSet)
                predictions = modelVisitor.predict(testSet.reactions.all(), whitelist, verbose=verbose)
                if verbose:
                    print "\t...finished predicting. Storing predictions...",
                newResDict = self._storePredictionComponents(predictions, statsModel)
    
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
                    print "predictions stored."
                    for response in self.outcomeDescriptors:
                        predDesc = response.predictedDescriptorType.objects.get(modelContainer=self, statsModel=statsModel, predictionOf=response)
                        conf_mtrx = predDesc.getConfusionMatrix()
                        
                        print "Confusion matrix for {}:".format(predDesc.heading)
                        print confusionMatrixString(conf_mtrx)
                        print "Accuracy: {:.3}".format(accuracy(conf_mtrx))
                        print "BCR: {:.3}".format(BCR(conf_mtrx))
            
            elif verbose:
                print "Test set is empty."

            if verbose:
                num_finished += 1
                end_time = datetime.datetime.now()
                elapsed = (end_time - overall_start_time)
                expected_finish = datetime.timedelta(seconds=(elapsed.total_seconds()*(num_models/float(num_finished)))) + overall_start_time
                print "{}. {} of {} models built.".format(end_time, num_finished, num_models)
                print "Elapsed model building time: {}. Expected completion time: {}".format(elapsed, expected_finish)

        if resDict:
            if verbose:
                print "Storing overall model predictions...",
            self._storePredictions(resDict)
            if verbose:
                print "Predictions stored"

        self.built = True
        if verbose:
            overall_end_time = datetime.datetime.now()
            print "Finished at {}".format(overall_end_time)

    def _storePredictionComponents(self, predictions, statsModel, resDict=None):
        """
        returns resDict, a dictionary of dictionaries of dictionaries
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
                vals.append(val)
            if vals:
                type(vals[0]).objects.bulk_create(vals)
        return resDict

    def _storePredictions(self, resDict):
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
                    #estimate the weighted average of the estimates from the component models
                    values = tuple(value for value in resDict[reaction][response].keys())
                    weights = tuple(weight for value, weight in resDict[reaction][response].items())
                    val = predDesc.createValue(reaction, average(values, weights=weights))
                    num_vals.append(val)
                else:
                    #find the 'competitors' with the highest number of votes and pick one at random (in case multiple categories have equal numbers of votes)
                    maxVotes = max(count for response, count in outcomeDict.items())
                    winners = [(response, count) for response, count in outcomeDict.items() if count == maxVotes]
                    winner = random.choice(winners)
                    val = predDesc.createValue(reaction, winner[0])
                    if isinstance(val, BoolRxnDescriptorValue):
                        bool_vals.append(val)
                    elif isinstance(val, OrdRxnDescriptorValue):
                        ord_vals.append(val)
                    elif isinstance(val, CatRxnDescriptorValue):
                        cat_vals.append(val)
                    else:
                        raise ValueError("Value just created is of unexpected type {}".format(type(val)))

                if response not in finalPredictions:
                    finalPredictions[response] = []
                finalPredictions[response].append( (reaction, val) )
                
        BoolRxnDescriptorValue.objects.bulk_create(bool_vals)
        NumRxnDescriptorValue.objects.bulk_create(num_vals)
        OrdRxnDescriptorValue.objects.bulk_create(ord_vals)
        CatRxnDescriptorValue.objects.bulk_create(cat_vals)
        return finalPredictions

    def predict(self, reactions):
        if self.built:
            resDict = {}

            for model in self.statsmodel_set:
                modelVisitor = getattr(visitorModules[self.modelVisitorLibrary], self.modelVisitorTool)(model)
                predictions = modelVisitor.predict(reactions)
                newResDict = self._storePredictionComponents(predictions, model)

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

            return self._storePredictions(resDict)
        else:
            raise RuntimeError('A model container cannot be used to make predictions before the build method has been called')

    @transaction.atomic
    def save(self, *args, **kwargs):
        super(ModelContainer, self).save(*args, **kwargs)

    def getOverallConfusionMatrices(self):
        confusion_matrix_list = []
        for descriptor in self.predictsDescriptors:
            if descriptor.statsModel is None:
                confusion_matrix_list.append( (descriptor.csvHeader, descriptor.getConfusionMatrix()) )
        return confusion_matrix_list


    def getConfusionMatrices(self):
        """
        Returns a list of lists of tuples of confusion matrices.
        Each entry of the outer list is for a different model.
        The first is the overall model, the rest component statsModels.
        For each model there is a list of tuples.
        Each tuple is of the form (descriptor_heading, confusion matrix)
        """

        confusion_matrix_lol = [self.getOverallConfusionMatrices()]

        # Retrieve the matrix for each component.
        for model in self.statsmodel_set.all():
            confusion_matrix_list = []
            for descriptor in self.predictsDescriptors:
                if descriptor.statsModel == model:
                    confusion_matrix_list.append( (descriptor.csvHeader, descriptor.getConfusionMatrix()) )
            confusion_matrix_lol.append(confusion_matrix_list)
            
        return confusion_matrix_lol

