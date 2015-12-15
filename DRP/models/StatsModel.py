"""A module containing only the StatsModel class."""
from django.db import models
from StatsModelTag import StatsModelTag
from ModelContainer import ModelContainer
from rxnDescriptors import NumRxnDescriptor, OrdRxnDescriptor, BoolRxnDescriptor, CatRxnDescriptor


class StatsModel(models.Model):

    """A class for describing a statistical model."""

    class Meta:
        app_label = 'DRP'

    fileName = models.FileField(upload_to='models', max_length=200, blank=True)
    """The filename in which this model is stored"""
    startTime = models.DateTimeField(default=None, null=True)
    endTime = models.DateTimeField(default=None, null=True)
    trainingSet = models.ForeignKey(DataSet)
    testSets = models.ManyToManyField(DataSet)
    container = models.ForeignKey(ModelContainer)

    # these fields are for use if a model should become invalidated
    invalid = models.BooleanField(default=False)
    regenerationOf = models.ForeignKey("self", blank=True, null=True, default=None)

    def invalidate(self):
        """Invalidate and regenerate the model instance."""
        self.invalid = True
        # TODO: populate other parts of this method
        self.save()
