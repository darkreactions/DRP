'''A module containing only the StatsModelTag class'''
from django.db import models

class StatsModelTag(models.Model):
  '''A class for tagging statsmodels with short text tags'''

  class Meta:
    app_label='DRP'

  text=models.CharField(max_length=200, unique=True)
