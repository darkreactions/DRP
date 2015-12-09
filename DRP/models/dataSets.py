'''A module containing classes to provide deletion protection to Performed reactions
used in StatsModels'''

from django.db import models
from PerformedReaction import PerformedReaction
from StatsModel import StatsModel

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
