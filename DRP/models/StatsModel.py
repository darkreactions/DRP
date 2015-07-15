from django.db import models
from Descriptor import Descriptor
from DRP.settings import LIBRARY_CHOICES, TOOL_CHOICES

class StatsModel(models.Model):

  class Meta:
    app_label='DRP'

  responses=models.ManyToManyField(DescriptorValue)
  dataSets=models.ManyToManyField(ConcreteReaction, through='DataSet')
  fileName=models.CharField('File for running and retrieval of model', max_length=200)
  description=models.TextField()
  active=models.BooleanField('Is this the active model?')
  library=models.CharField(
  tool=models.CharField(
  start_time=models.DateTimeField()
  end_time=models.DataTimeField(default=None, null=True)
  iterations=models.IntegerField()
  library=models.CharField(max_length=200, choices=LIBRARY_CHOICES)
  tool=models.Charfield(max_length=200, choices=TOOL_CHOICES)
