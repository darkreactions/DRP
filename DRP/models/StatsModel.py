'''A module containing only the StatsModel class'''
from django.db import models
from StatsModelTag import StatsModelTag
from descriptors import Descriptor 
from django.conf import settings

class StatsModel(models.Model):
  '''A class for describing a statistical model generated to describe statistical models generated as
  a part of the DRP'''

  class Meta:
    app_label='DRP'

  fileName=models.CharField('File for running and retrieval of model', max_length=200)
  '''The filename in which this model is stored'''
  description=models.TextField()
  active=models.BooleanField('Is this the active model?')
  start_time=models.DateTimeField()
  end_time=models.DateTimeField(default=None, null=True)
  iterations=models.IntegerField()
  library=models.CharField(max_length=200, choices=settings.LIBRARY_CHOICES)
  tool=models.CharField(max_length=200, choices=settings.TOOL_CHOICES)
  tags=models.ManyToManyField(StatsModelTag)
  descriptors=models.ManyToManyField(Descriptor) 
