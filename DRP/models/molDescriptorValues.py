'''A module containign only the DescriptorValue class'''
from django.db import models
from descriptorValues import CategoricalDescriptorValue, BooleanDescriptorValue, NumericDescriptorValue, OrdinalDescriptorValue
#from Compound import DRP.Compound - retain this line for clarity
from django.core.exceptions import ValidationError

class CatMolDescriptorValue(CategoricalDescriptorValue):
  '''Contains the value of a categorical descriptor for a compound'''

  class Meta:
    app_label="DRP"
    verbose_name='Categorical Molecular Descriptor Value'
    unique_together=('descriptor', 'compound')

  compound = models.ForeignKey('DRP.Compound')

class BoolMolDescriptorValue(BooleanDescriptorValue):
  '''Contains the value of a boolean descriptor for a compound'''

  class Meta:
    app_label="DRP"
    verbose_name='Boolean Molecular Descriptor Value'
    unique_together=('descriptor', 'compound')

  compound = models.ForeignKey('DRP.Compound')

class NumMolDescriptorValue(NumericDescriptorValue):
  '''Contains the numeric value of a descriptor for a compound'''

  class Meta:
    app_label="DRP"
    verbose_name='Numeric Molecular Descriptor Value'
    unique_together=('descriptor', 'compound')

  compound = models.ForeignKey('DRP.Compound')

class OrdMolDescriptorValue(OrdinalDescriptorValue):
  '''Contains the ordinal value of a descriptor for a compound'''

  class Meta:
    app_label="DRP"
    verbose_name='Ordinal Molecular Descriptor Value'
    unique_together=('descriptor', 'compound')

  compound = models.ForeignKey('DRP.Compound')
