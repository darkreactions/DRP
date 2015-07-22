'''A module containing only the StatsModelTag class'''
from django.db import models

class StatsModelTags(models.Model)
  '''A class for tagging statsmodels with short text tags'''

  class Meta:
    app_label='DRP'

  text=models.Charfield(max_length=200, unique=True)
