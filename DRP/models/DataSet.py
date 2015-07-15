from django.db import models
from StatsModel import StatsModel

class DataSet(models.Model):

  Meta:
    app_label='DRP'

  reaction=models.ForeignKey(Reaction)
  model=models.ForeignKey(StatsModel)
  isTestSet=models.BooleanField()
  isTrainingSet=models.BooleanField()
