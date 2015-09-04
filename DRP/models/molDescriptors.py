'''A module containing Classes permitting the representation of molecular descriptors'''
from django.db import models
from django.core.exceptions import ValidationError
from descriptors import Descriptor, CategoricalDescriptor, OrdinalDescriptor, BooleanDescriptor
from descriptors import CategoricalDescriptorPermittedValue, NumericDescriptor
from django.core.validators import RegexValidator


class CatMolDescriptor(CategoricalDescriptor):
  '''A class which describes a categorical molecular descriptors'''

  class Meta:
    app_label='DRP'
    verbose_name= 'Categorical Molecular Descriptor'

class OrdMolDescriptor(OrdinalDescriptor):
  '''A class which represents an ordinal descriptor'''
  
  class Meta:
    verbose_name= 'Ordinal Molecular Descriptor'
    app_label='DRP'

class NumMolDescriptor(NumericDescriptor):
  '''A class which represents a numerical descriptor'''

  class Meta:
    app_label='DRP'
    verbose_name= 'Numerical Molecular Descriptor'


class BoolMolDescriptor(BooleanDescriptor):
  '''A class which represents a boolean descriptors'''

  class Meta:
    app_label='DRP'
    verbose_name= 'Boolean Molecular Descriptor'

