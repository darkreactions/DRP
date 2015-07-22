'''A module containing only the RecommendedReaction class'''
from django.db import models
from Reaction import Reaction
import StatsModel

class RecommendedReaction(Reaction):

  class Meta:
    app_label='DRP'

  score=models.FloatField()
  model=models.ForeignKey(StatsModel.StatsModel, null=True)
  seed=models.ForeignKey(Reaction, null=True)
  nonsense=models.BooleanField()
  hidden=models.BooleanField()
  saved=models.BooleanField()
  reference=models.CharField('Text Reference', max_length=200)
