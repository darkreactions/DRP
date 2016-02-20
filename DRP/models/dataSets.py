'''A module containing classes to provide deletion protection to Performed reactions
used in StatsModels. Whilst the DataSet model may, to the uninitiated, appear to be
an erroneous proxy for a many-to-many relationship between Reactions and Models,
This allows the datasets to exist independently of the models.'''

from django.db import models
from PerformedReaction import PerformedReaction


class DataSet(models.Model):

  class Meta:
    app_label="DRP"

  name = models.CharField(max_length=200, unique=True)
  reactions = models.ManyToManyField(PerformedReaction, through="DataSetRelation")


class DataSetRelation(models.Model):

  class Meta:
    app_label="DRP"
    unique_together = ("dataSet", "reaction")

  reaction = models.ForeignKey(PerformedReaction, on_delete=models.PROTECT)
  dataSet = models.ForeignKey(DataSet)
