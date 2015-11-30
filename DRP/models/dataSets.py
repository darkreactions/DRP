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

  name = models.CharField(max_length=200, unique=True)
  model = models.ForeignKey(StatsModel)
  reactions = models.ManyToManyField(PerformedReaction, through="TestSetRelation")


class TestSetRelation(models.Model):

  class Meta:
    app_label="DRP"
    unique_together = ("test_set", "reaction")

  reaction = models.ForeignKey(PerformedReaction, on_delete=models.PROTECT)
  test_set = models.ForeignKey(TestSet)



