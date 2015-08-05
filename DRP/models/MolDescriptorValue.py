'''A module containign only the DescriptorValue class'''
from django.db import models
from MolDescriptor import MolDescriptor
from Reaction import Reaction
from Compound import Compound
from StatsModel import StatsModel

class MolDescriptorValue(models.Model):
  '''Contains Relationships between Compounds and their descriptors'''

  class Meta:
    app_label="DRP"
    verbose_name='Molecular Descriptor Value'
    unique_together=('descriptor', 'compound')

  descriptor = models.ForeignKey(MolDescriptor)
  compound = models.ForeignKey(Compound, null=True, unique=False, default=None)
  booleanValue= models.NullBooleanField('Value if descriptor is a boolean', null=True)
  ordValue = models.PositiveIntegerField('Value if descriptor is an ordinal', null=True)
  catValue = models.CharField('Value if descriptor is a category', max_length=200, null=True)
  numValue = models.FloatField('Value if descritpor is continuous', null=True)
