"""A module containing only the StatsModel class."""
from django.db import models
from StatsModelTag import StatsModelTag
from ModelContainer import ModelContainer
from rxnDescriptors import NumRxnDescriptor, OrdRxnDescriptor, BoolRxnDescriptor, CatRxnDescriptor

class DescriptorAttribute(object):

    def __get__(self, model, modelType=None):
        return chain(model.boolRxnDescriptors.all(), model.ordRxnDescriptors.all(), model.catRxnDescriptors.all(), model.numRxnDescriptors.all())

    def __set__(self, model, descriptors):
        model.boolRxnDescriptors.clear()
        model.ordRxnDescriptors.clear()
        model.catRxnDescriptors.clear()
        model.numRxnDescriptors.clear()
        for descriptor in descriptors:
            if isinstance(descriptor, BoolRxnDescriptor):
                model.boolRxnDescriptors.add(descriptor):
            elif isinstance(descriptor, OrdRxnDescriptor):
                model.ordRxnDescriptors.add(descriptor)
            elif isinstance(descriptor, CatRxnDescriptor):
                model.catRxnDescriptors.add(CatRxnDescriptor)
            elif isinstance(descriptor, NumRxnDescriptor):
                model.numRxnDescriptors.add(NumRxnDescriptor)
            else:
                raise ValueError('An invalid object was assigned as a descriptor')

    def __delete__(self, model):
        model.boolRxnDescriptors.clear()
        model.numRxnDescriptors.clear()
        model.catRxnDescriptors.clear()
        model.ordRxnDescriptors.clear()

class OutComeDescriptorAttribute(object):

    def __get__(self, model, modelType=None):
        return chain(model.outcomeBoolRxnDescriptors.all(), model.outcomeOrdRxnDescriptors.all(), model.outcomeCatRxnDescriptors.all(), model.outcomeNumRxnDescriptors.all())

    def __set__(self, model, descriptors):
        model.outcomeBoolRxnDescriptors.clear()
        model.outcomeOrdRxnDescriptors.clear()
        model.outcomeCatRxnDescriptors.clear()
        model.outcomeNumRxnDescriptors.clear()
        for descriptor in descriptors:
            if isinstance(descriptor, BoolRxnDescriptor):
                model.outcomeBoolRxnDescriptors.add(descriptor):
            elif isinstance(descriptor, OrdRxnDescriptor):
                model.outcomeOrdRxnDescriptors.add(descriptor)
            elif isinstance(descriptor, CatRxnDescriptor):
                model.outcomeCatRxnDescriptors.add(CatRxnDescriptor)
            elif isinstance(descriptor, NumRxnDescriptor):
                model.outcomeNumRxnDescriptors.add(NumRxnDescriptor)
            else:
                raise ValueError('An invalid object was assigned as a descriptor')

    def __delete__(self, model):
        model.outcomeBoolRxnDescriptors.clear()
        model.outcomeNumRxnDescriptors.clear()
        model.outcomeCatRxnDescriptors.clear()
        model.outcomeOrdRxnDescriptors.clear()

class StatsModel(models.Model):

    """A class for describing a statistical model."""

    class Meta:
        app_label = 'DRP'

    fileName = models.FileField(upload_to='models', max_length=200)
    """The filename in which this model is stored"""
    description = models.TextField()
    active = models.BooleanField('Is this the active model?', default=False)
    start_time = models.DateTimeField(default=None, null=True)
    end_time = models.DateTimeField(default=None, null=True)
    iterations = models.IntegerField()
    tags = models.ManyToManyField(StatsModelTag)

    container = models.ForeignKey(ModelContainer)

    descriptors = DescriptorAttribute()
    boolRxnDescriptors = models.ManyToManyField(BoolRxnDescriptor)
    ordRxnDescriptors = models.ManyToManyField(OrdRxnDescriptor)
    catRxnDescriptors = models.ManyToManyField(CatRxnDescriptor)
    numRxnDescriptors = models.ManyToManyField(NumRxnDescriptor)
    """The input descriptors for the model."""

    outcomeDescriptors OutComeDescriptorAttribute()
    outcomeBoolRxnDescriptors = models.ManyToManyField(BoolRxnDescriptor, related_name='outcomeForModels'))
    outcomeOrdRxnDescriptors = models.ManyToManyField(OrdRxnDescriptor, related_name='outcomeForModels'))
    outcomeCatRxnDescriptors = models.ManyToManyField(CatRxnDescriptor, related_name='outcomeForModels'))
    outcomeNumRxnDescriptors = models.ManyToManyField(NumRxnDescriptor, related_name='outcomeForModels'))
    """The descriptors which are being used as outcomes for this model.

    For models which make predictions about descriptors, it is probably
    most appropriate to make the descriptor label for the predicted
    descriptor related to the model, for instance for an outcome
    descriptor called "outcome", you might consider
    "outcome_predicted_by_model_id_1", where 1 is the model primary
    key.
    """

    predictsDescriptors = models.ManyToManyField(
        Descriptor, related_name="predictedByModels")
    """The descriptors which this model predicts values for."""

    # these fields are for use if a model should become invalidated
    invalid = models.BooleanField(default=False)
    regenerationOf = models.ForeignKey("self", blank=True, null=True, default=None)
    snapShot = models.FileField(
        upload_to='model_snapshots',
        max_length=200,
        null=True,
        default=None)

    def invalidate(self):
        """Invalidate and regenerate the model instance."""
        self.invalid = True
        # TODO: populate other parts of this method
        self.save()
