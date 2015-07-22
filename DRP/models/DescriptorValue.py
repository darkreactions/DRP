'''A module containign only the DescriptorValue class'''
from django.db import models
from Descriptor import Descriptor
from Reaction import Reaction
from Compound import Compound
from StatsModel import StatsModel

class DescriptorValue(models.Model):
  '''Contains Relationships between both Reactions and Compounds and their descriptors'''

  class Meta:
    app_label="DRP"

  descriptor = models.ForeignKey(Descriptor)
  Reaction = models.ForeignKey(Reaction, null=True, unique=False, default=None)
  Compound = models.ForeignKey(Compound, null=True, unique=False, default=None)
  booleanValue= models.BooleanField('Value if descriptor is a boolean', null=True)
  ordValue = models.PositiveIntegerField('Value if descriptor is an ordinal', null=True)
  catValue = models.CharField('Value if descriptor is a category', max_length=200, null=True)
  numValue = models.FloatField('Value if descritpor is continuous', null=True)
  isPredicted=models.BooleanField()
  '''Some stored descriptors will be predicted outcomes formed from the basis of a statistical model.
  This provides that indication
  '''
  model=models.ForeignKey(StatsModel, unique=False, null=True, default=None)
  '''If this value was predicted by a statistical model, reference that model'''
