'''A module containign only the DescriptorValue class'''
from django.db import models
from Reaction import Reaction
from descriptorValues import CategoricalDescriptorValue, OrdinalDescriptorValue,BooleanDescriptorValue, NumericDescriptorValue  
from StatsModel import StatsModel

class RxnDescriptorValue(models.Model):
  '''Contains Relationships between Reactions and their descriptors'''

  class Meta:
    app_label="DRP"
    abstract=True

  reaction = models.ForeignKey(Reaction, unique=False)
  model=models.ForeignKey(StatsModel, unique=False, null=True, default=None)
  '''If this value was predicted by a statistical model, reference that model'''

class CatRxnDescriptorValue(CategoricalDescriptorValue, RxnDescriptorValue):
  '''Contains the value of a categorical descriptor for a reaction'''

  class Meta:
    app_label="DRP"
    verbose_name='Categorical Reaction Descriptor Value'
    unique_together=('reaction', 'descriptor')

class BoolRxnDescriptorValue(BooleanDescriptorValue, RxnDescriptorValue):
  '''Contains the value of a boolean descriptor for a reaction'''

  class Meta:
    app_label="DRP"
    verbose_name='Boolean Reaction Descriptor Value'
    unique_together=('reaction', 'descriptor')

class NumRxnDescriptorValue(NumericDescriptorValue, RxnDescriptorValue):
  '''Contains the numeric value of a descriptor for a reaction'''

  class Meta:
    app_label="DRP"
    verbose_name='Numeric Reaction Descriptor Value'
    unique_together=('reaction', 'descriptor')

class OrdRxnDescriptorValue(OrdinalDescriptorValue, RxnDescriptorValue):
  '''Contains the ordinal value of a descriptor for a reaction'''

  class Meta:
    app_label="DRP"
    verbose_name='Ordinal Reaction Descriptor Value'
    unique_together=('reaction', 'descriptor')
