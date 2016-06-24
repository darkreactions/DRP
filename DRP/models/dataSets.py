"""
A module containing classes to group reactions into sets.

These classses provide  deletion protection to Performed reactions
used in StatsModels. Whilst the DataSet model may, to the uninitiated, appear to be
an erroneous proxy for a many-to-many relationship between Reactions and Models,
This allows the datasets to exist independently of the models.
"""

from django.db import models
from .performedReaction import PerformedReaction


class DataSet(models.Model):

    """A set of reactions."""

    class Meta:
        app_label = "DRP"

    name = models.CharField(max_length=200, unique=True)
    reactions = models.ManyToManyField(
        PerformedReaction, through="DataSetRelation")

    # TODO: This belongs on a manager to be djangonic.
    @classmethod
    def create(cls, name, data):
        """Bulk create a set of datasetrelations."""
        dataSet = cls(name=name)
        dataSet.save()
        dsrs = [DataSetRelation(dataSet=dataSet, reaction=datum)
                for datum in data]
        DataSetRelation.objects.bulk_create(dsrs)

        return dataSet


class DataSetRelation(models.Model):

    """Defines the relationships between a data set and a reaction."""

    class Meta:
        app_label = "DRP"
        unique_together = ("dataSet", "reaction")

    reaction = models.ForeignKey(PerformedReaction, on_delete=models.PROTECT)
    dataSet = models.ForeignKey(DataSet)
