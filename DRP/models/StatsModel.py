"""A module containing only the StatsModel class."""
from dataSets import DataSet
from django.db import models


class StatsModel(models.Model):

    """A class for describing a statistical model."""

    class Meta:
        app_label = 'DRP'

    outputFile = models.FileField(upload_to='models', max_length=200, blank=True)
    inputFile = models.FileField(upload_to='model_inputs', max_length=255, blank=True)
    """The filename in which this model is stored"""
    startTime = models.DateTimeField(default=None, null=True)
    endTime = models.DateTimeField(default=None, null=True)
    trainingSet = models.ForeignKey(DataSet, related_name='trainingSetFor')
    testSets = models.ManyToManyField(DataSet, related_name='testSetsFor')
    container = models.ForeignKey("DRP.ModelContainer")

    # these fields are for use if a model should become invalidated
    invalid = models.BooleanField(default=False)
    regenerationOf = models.ForeignKey("self", blank=True, null=True, default=None)

    def invalidate(self):
        """Invalidate the model instance."""
        self.invalid = True
        self.save()
