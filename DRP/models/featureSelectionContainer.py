"""A module containing FeatureSelectionContainer class and related classes for descriptor attributes."""
from django.db import models
from django.conf import settings
from DRP.models import DataSet, Descriptor, BoolRxnDescriptor, OrdRxnDescriptor, NumRxnDescriptor, CatRxnDescriptor
import importlib
from itertools import chain
import datetime
import json
import logging

logger = logging.getLogger(__name__)

featureVisitorModules = {library: importlib.import_module(
    settings.FEATURE_SELECTION_LIBS_DIR + "." + library) for library in settings.FEATURE_SELECTION_LIBS}

FEATURE_SELECTION_TOOL_CHOICES = tuple(
    tool for library in featureVisitorModules.values() for tool in library.tools)


class DescriptorAttribute(object):
    """An attribute manager object which allows the setting and deletion of the related descriptors."""

    def __get__(self, featureSelectionContainer, featureSelectionContainerType=None):
        """Return a chain of each of the querysets."""
        return chain(featureSelectionContainer.boolRxnDescriptors.all(), featureSelectionContainer.ordRxnDescriptors.all(), featureSelectionContainer.catRxnDescriptors.all(), featureSelectionContainer.numRxnDescriptors.all())

    def __set__(self, featureSelectionContainer, descriptors):
        """Clear the descriptor sets and then set them equal to the provided queryset."""
        featureSelectionContainer.boolRxnDescriptors.clear()
        featureSelectionContainer.ordRxnDescriptors.clear()
        featureSelectionContainer.catRxnDescriptors.clear()
        featureSelectionContainer.numRxnDescriptors.clear()
        for descriptor in descriptors:
            desc = None
            try:
                desc = BoolRxnDescriptor.objects.get(id=descriptor.id)
                featureSelectionContainer.boolRxnDescriptors.add(desc)
            except BoolRxnDescriptor.DoesNotExist:
                pass

            try:
                desc = OrdRxnDescriptor.objects.get(id=descriptor.id)
                featureSelectionContainer.ordRxnDescriptors.add(desc)
            except OrdRxnDescriptor.DoesNotExist:
                pass

            try:
                desc = CatRxnDescriptor.objects.get(id=descriptor.id)
                featureSelectionContainer.catRxnDescriptors.add(desc)
            except CatRxnDescriptor.DoesNotExist:
                pass

            try:
                desc = NumRxnDescriptor.objects.get(id=descriptor.id)
                featureSelectionContainer.numRxnDescriptors.add(desc)
            except NumRxnDescriptor.DoesNotExist:
                pass

            if desc is None:
                raise ValueError(
                    'An invalid object with heading {} and type {} was assigned as a descriptor'.format(descriptor.heading, type(descriptor)))

    def __delete__(self, featureSelectionContainer):
        """Simply clear the querysets."""
        featureSelectionContainer.boolRxnDescriptors.clear()
        featureSelectionContainer.numRxnDescriptors.clear()
        featureSelectionContainer.catRxnDescriptors.clear()
        featureSelectionContainer.ordRxnDescriptors.clear()


class ChosenDescriptorAttribute(object):
    """An attribute manager object which allows the setting and deletion of the related descriptors chosen by feature selection."""

    def __get__(self, featureSelectionContainer, featureSelectionContainerType=None):
        """Return a chain of each of the querysets."""
        return chain(featureSelectionContainer.chosenBoolRxnDescriptors.all(), featureSelectionContainer.chosenOrdRxnDescriptors.all(), featureSelectionContainer.chosenCatRxnDescriptors.all(), featureSelectionContainer.chosenNumRxnDescriptors.all())

    def __set__(self, featureSelectionContainer, descriptors):
        """Clear the descriptor sets and then set them equal to the provided queryset."""
        featureSelectionContainer.chosenBoolRxnDescriptors.clear()
        featureSelectionContainer.chosenOrdRxnDescriptors.clear()
        featureSelectionContainer.chosenCatRxnDescriptors.clear()
        featureSelectionContainer.chosenNumRxnDescriptors.clear()
        for descriptor in descriptors:
            desc = None
            try:
                desc = BoolRxnDescriptor.objects.get(id=descriptor.id)
                featureSelectionContainer.chosenBoolRxnDescriptors.add(desc)
            except BoolRxnDescriptor.DoesNotExist:
                pass

            try:
                desc = OrdRxnDescriptor.objects.get(id=descriptor.id)
                featureSelectionContainer.chosenOrdRxnDescriptors.add(desc)
            except OrdRxnDescriptor.DoesNotExist:
                pass

            try:
                desc = CatRxnDescriptor.objects.get(id=descriptor.id)
                featureSelectionContainer.chosenCatRxnDescriptors.add(desc)
            except CatRxnDescriptor.DoesNotExist:
                pass

            try:
                desc = NumRxnDescriptor.objects.get(id=descriptor.id)
                featureSelectionContainer.chosenNumRxnDescriptors.add(desc)
            except NumRxnDescriptor.DoesNotExist:
                pass

            if desc is None:
                raise ValueError(
                    'An invalid object with heading {} and type {} was assigned as a descriptor'.format(descriptor.heading, type(descriptor)))

    def __delete__(self, featureSelectionContainer):
        """Simply clear the querysets."""
        featureSelectionContainer.chosenBoolRxnDescriptors.clear()
        featureSelectionContainer.chosenNumRxnDescriptors.clear()
        featureSelectionContainer.chosenCatRxnDescriptors.clear()
        featureSelectionContainer.chosenOrdRxnDescriptors.clear()


class OutcomeDescriptorAttribute(object):
    """An attribute manager object which allows the setting and deletion of the related outcome descriptors."""

    def __get__(self, featureSelectionContainer, featureSelectionContainerType=None):
        """Return a chain of each of the querysets."""
        return chain(featureSelectionContainer.outcomeBoolRxnDescriptors.all(), featureSelectionContainer.outcomeOrdRxnDescriptors.all(), featureSelectionContainer.outcomeCatRxnDescriptors.all(), featureSelectionContainer.outcomeNumRxnDescriptors.all())

    def __set__(self, featureSelectionContainer, descriptors):
        """Clear the descriptor sets and then set them equal to the provided queryset."""
        featureSelectionContainer.outcomeBoolRxnDescriptors.clear()
        featureSelectionContainer.outcomeOrdRxnDescriptors.clear()
        featureSelectionContainer.outcomeCatRxnDescriptors.clear()
        featureSelectionContainer.outcomeNumRxnDescriptors.clear()
        for descriptor in descriptors:
            desc = None
            try:
                desc = BoolRxnDescriptor.objects.get(id=descriptor.id)
                featureSelectionContainer.outcomeBoolRxnDescriptors.add(desc)
            except BoolRxnDescriptor.DoesNotExist:
                pass

            try:
                desc = OrdRxnDescriptor.objects.get(id=descriptor.id)
                featureSelectionContainer.outcomeOrdRxnDescriptors.add(desc)
            except OrdRxnDescriptor.DoesNotExist:
                pass

            try:
                desc = CatRxnDescriptor.objects.get(id=descriptor.id)
                featureSelectionContainer.outcomeCatRxnDescriptors.add(desc)
            except CatRxnDescriptor.DoesNotExist:
                pass

            try:
                desc = NumRxnDescriptor.objects.get(id=descriptor.id)
                featureSelectionContainer.outcomeNumRxnDescriptors.add(desc)
            except NumRxnDescriptor.DoesNotExist:
                pass

            if desc is None:
                raise ValueError(
                    'An invalid object was assigned as a descriptor')

    def __delete__(self, featureSelectionContainer):
        """Simply clear the querysets."""
        featureSelectionContainer.outcomeBoolRxnDescriptors.clear()
        featureSelectionContainer.outcomeNumRxnDescriptors.clear()
        featureSelectionContainer.outcomeCatRxnDescriptors.clear()
        featureSelectionContainer.outcomeOrdRxnDescriptors.clear()


class FeatureSelectionContainer(models.Model):
    """A class for encapsulating a feature selection model."""

    description = models.TextField(default='', blank=True)
    featureVisitorLibrary = models.CharField(
        max_length=200, default='', blank=True)
    featureVisitorTool = models.CharField(
        max_length=200, default='', blank=True)
    featureVisitorOptions = models.TextField(
        null=False, blank=True, default="{}")
    startTime = models.DateTimeField(default=None, null=True, blank=True)
    endTime = models.DateTimeField(default=None, null=True, blank=True)
    trainingSet = models.ForeignKey(
        DataSet, related_name='trainingSetForFeatureSelection', null=True)
    built = models.BooleanField(
        'Has the build procedure been called with this container?', editable=False, default=False)

    descriptors = DescriptorAttribute()
    boolRxnDescriptors = models.ManyToManyField(BoolRxnDescriptor)
    ordRxnDescriptors = models.ManyToManyField(OrdRxnDescriptor)
    catRxnDescriptors = models.ManyToManyField(CatRxnDescriptor)
    numRxnDescriptors = models.ManyToManyField(NumRxnDescriptor)

    outcomeDescriptors = OutcomeDescriptorAttribute()
    outcomeBoolRxnDescriptors = models.ManyToManyField(
        BoolRxnDescriptor, related_name='outcomeForFeatureSelection')
    outcomeOrdRxnDescriptors = models.ManyToManyField(
        OrdRxnDescriptor, related_name='outcomeForFeatureSelections')
    outcomeCatRxnDescriptors = models.ManyToManyField(
        CatRxnDescriptor, related_name='outcomeForFeatureSelections')
    outcomeNumRxnDescriptors = models.ManyToManyField(
        NumRxnDescriptor, related_name='outcomeForFeatureSelections')

    chosenDescriptors = ChosenDescriptorAttribute()
    chosenBoolRxnDescriptors = models.ManyToManyField(
        BoolRxnDescriptor, related_name='chosenForFeatureSelection')
    chosenOrdRxnDescriptors = models.ManyToManyField(
        OrdRxnDescriptor, related_name='chosenForFeatureSelection')
    chosenCatRxnDescriptors = models.ManyToManyField(
        CatRxnDescriptor, related_name='chosenForFeatureSelection')
    chosenNumRxnDescriptors = models.ManyToManyField(
        NumRxnDescriptor, related_name='chosenForFeatureSelection')

    @classmethod
    def create(cls, featureVisitorLibrary, featureVisitorTool, predictors, responses, featureVisitorOptions=None, reactions=None, trainingSet=None, description=""):
        """
        Create a feature selection container with a set of options.

        featureVisitorLibrary sepcifies the library/plugin being used, while the featureVisitorTool should specify the actual tool.
        predictors is the descriptors being used as independent variables, responses is a descriptor queryset for the independent variables.
        reactions or training sets should be specified, and this will define how the training data is defined.
        """
        container = cls(featureVisitorLibrary=featureVisitorLibrary,
                        featureVisitorTool=featureVisitorTool, description=description)
        container.save()  # need pk

        if featureVisitorOptions is not None:
            container.featureVisitorOptions = json.dumps(featureVisitorOptions)

        if trainingSet is None:
            assert(reactions is not None)
            container.trainingSet = DataSet.create('{}_{}_{}'.format(
                container.featureVisitorLibrary, container.featureVisitorLibrary, container.pk), reactions)
        else:
            container.trainingSet = trainingSet

        container.descriptors = predictors
        container.outcomeDescriptors = responses

        return container

    def build(self, verbose=False):
        """Take all options confirmed so far and generates a feature set."""
        if self.built:
            raise RuntimeError(
                "Cannot build a feature selection that has already been built.")

        visitorOptions = json.loads(self.featureVisitorOptions)
        featureVisitor = getattr(featureVisitorModules[
                                 self.featureVisitorLibrary], self.featureVisitorTool)(container=self, **visitorOptions)

        self.startTime = datetime.datetime.now()

        if verbose:
            logger.info("{}, training on {} reactions...".format(self.startTime, self.trainingSet.reactions.count()))
        chosen_descriptor_headers = featureVisitor.train(verbose=verbose)

        self.endTime = datetime.datetime.now()
        if verbose:
            logging.info("\t...Trained. Finished at {}.".format(self.endTime))

        # TODO XXX Make this not hacky af
        chosen_descriptor_list = []
        for desc in self.descriptors:
            if desc.csvHeader in chosen_descriptor_headers:
                chosen_descriptor_list.append(desc)

        self.chosenDescriptors = chosen_descriptor_list
        self.built = True
