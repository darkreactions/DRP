"""A module containing only the StatsModel class."""
from django.db import models
from StatsModelTag import StatsModelTag
from ModelContainer import ModelContainer
from descriptors import Descriptor


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

    descriptors = models.ManyToManyField(Descriptor)
    #TODO: set these as Attributes as per hte model visitor classes, and have 'real' relationships set
    # up to each descriptor type
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
