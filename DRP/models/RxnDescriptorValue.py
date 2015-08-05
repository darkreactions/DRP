'''A module containign only the DescriptorValue class'''
from django.db import models
from RxnDescriptor import RxnDescriptor
from Reaction import Reaction
from Compound import Compound
from StatsModel import StatsModel

class RxnDescriptorValue(models.Model):
  '''Contains Relationships between Reactions and their descriptors'''

  class Meta:
    app_label="DRP"
    verbose_name='Reaction Descriptor Value'
    unique_together=('reaction', 'descriptor')

  descriptor = models.ForeignKey(MolDescriptor)
  reaction = models.ForeignKey(Reaction, null=True, unique=False, default=None)
  booleanValue= models.NullBooleanField('Value if descriptor is a boolean', null=True)
  ordValue = models.PositiveIntegerField('Value if descriptor is an ordinal', null=True)
  catValue = models.CharField('Value if descriptor is a category', max_length=200, null=True)
  numValue = models.FloatField('Value if descritpor is continuous', null=True)
  isPredicted=models.BooleanField()
  '''Some stored descriptors will be predicted outcomes formed from the basis of a statistical model.
  This provides that indication
  '''
  model=models.ForeignKey(StatsModel, unique=False, null=True, default=None)
  '''If this value was predicted by a statistical model, reference that model'''
