'''A module containing the reactions descriptors'''
from django.db import models
from descriptors import Descriptor, CategoricalDescriptor, OrdinalDescriptor, BooleanDescriptor
from descriptors import CategoricalDescriptorPermittedValue, NumericDescriptor

class CatRxnDescriptor(models.Model):
  '''A class which describes a descriptor- a value which describes a system such as a compound or a reaction'''
  
  class Meta:
    app_label='DRP'
    verbose_name = 'Categorical Reaction Descriptor'

class OrdRxnDescriptor(OrdinalDescriptor):
  '''A class which represents an ordinal descriptor'''
  
  class Meta:
    verbose_name= 'Ordinal Reaction Descriptor'
    app_label='DRP'

class NumRxnDescriptor(NumericDescriptor):
  '''A class which represents a numerical descriptor'''

  class Meta:
    app_label='DRP'
    verbose_name= 'Numerical Reaction Descriptor'


class BoolRxnDescriptor(BooleanDescriptor):
  '''A class which represents a boolean descriptors'''

  class Meta:
    app_label='DRP'
    verbose_name= 'Boolean Reaction Descriptor'


