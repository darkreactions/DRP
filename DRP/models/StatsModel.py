'''A module containing only the StatsModel class'''
from django.db import models
from os import path
from StatsModelTag import StatsModelTag
from descriptors import Descriptor 
from django.conf import settings

class StatsModel(models.Model):


    """A class for describing a statistical model."""

    class Meta:
        app_label='DRP'

    fileName = models.FileField(upload_to='models', max_length=200)
    '''The filename in which this model is stored'''
    description = models.TextField()
    active = models.BooleanField('Is this the active model?')
    start_time = models.DateTimeField()
    end_time = models.DateTimeField(default=None, null=True)
    iterations = models.IntegerField()
    library = models.CharField(max_length=200, choices=settings.LIBRARY_CHOICES)
    tool = models.CharField(max_length=200, choices=settings.TOOL_CHOICES)
    tags = models.ManyToManyField(StatsModelTag)

    descriptors = models.ManyToManyField(Descriptor) 
    """The input descriptors for the model."""

    outcomeDescriptors = models.ManyToManyField(
        Descriptor, related_name='outcomeForModels')
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

    #these fields are for use if a model should become invalidated
    invalid = models.BooleanField(default=False)
    regenerationOf = models.ForeignKey(
            "self",
            blank=True,
            null=True,
            default=None)
    snapShot = models.FileField(
        upload_to='model_snapshots',
        max_length=200,
        null=True,
        default=None)

    def invalidate(self):
        self.invalid=True
        #TODO: populate other parts of this method
        self.save()
