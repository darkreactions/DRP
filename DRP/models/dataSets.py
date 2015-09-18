'''A module containing classes to provide deletion protection to Performed reactions
used in StatsModels'''

from django.db import models
from PerformedReaction import PerformedReaction
from StatsModel import StatsModel

class TrainingSet(models.Model):

  class Meta:
    app_label="DRP"

  reaction = models.ForeignKey(PerformedReaction, on_delete=models.PROTECT) 
  model = models.ForeignKey(StatsModel)

class TestSet(models.Model):

  class Meta:
    app_label="DRP"

  reaction = models.ForeignKey(PerformedReaction, on_delete=models.PROTECT) 
  model = models.ForeignKey(StatsModel)
