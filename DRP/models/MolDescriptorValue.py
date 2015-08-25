'''A module containign only the DescriptorValue class'''
from django.db import models
from MolDescriptor import CatMolDescriptor, BoolMolDescriptor, NumMolDescriptor, OrdMolDescriptor,CatMolDescriptorPermitted
from Reaction import Reaction
from Compound import Compound
from StatsModel import StatsModel

#TODO: implement methods which permit these to be used in the same way as their actual value attributes
#I.E. NumMolDescriptors can be multiplied and divided etc, returning the appropriate type.

#TODO: implment clean methods`

class CatMolDescriptorValue(models.Model):
  '''Contains the value of a categorical descriptor for a compound'''

  class Meta:
    app_label="DRP"
    verbose_name='Categorical Molecular Descriptor Value'
    unique_together=('descriptor', 'compound')

  descriptor = models.ForeignKey(CatMolDescriptor)
  compound = models.ForeignKey(Compound)
  value = models.ForeignKey(CatMolDescriptorPermitted)

  def __eq__(self, other):
    if isinstance(other, CatMolDescriptorValue):
      return self.value.value == other.value.value
    elif isinstance(other, CatMolDescriptorPermitted):
      return self.value.value == other.value
    else:
      return self.value.value == other

class BoolMolDescriptorValue(models.Model):
  '''Contains the value of a boolean descriptor for a compound'''

  class Meta:
    app_label="DRP"
    verbose_name='Boolean Molecular Descriptor Value'
    unique_together=('descriptor', 'compound')

  descriptor = models.ForeignKey(BoolMolDescriptor)
  compound = models.ForeignKey(Compound)
  value= models.NullBooleanField('Value if descriptor is a boolean', null=True)

  def __nonzero__(self):
    return self.value

class NumMolDescriptorValue(models.Model):
  '''Contains the numeric value of a descriptor for a compound'''

  class Meta:
    app_label="DRP"
    verbose_name='Boolean Molecular Descriptor Value'
    unique_together=('descriptor', 'compound')

  descriptor = models.ForeignKey(BoolMolDescriptor)
  compound = models.ForeignKey(Compound)
  value=models.FloatField(null=True)

  def __nonzero__(self):
    return bool(self.value)

  def __eq__(self, other):
    if isinstance(other, NumMolDescriptorValue):
      return self.value == other.value 
    else:
      return self.value == other

  def __gt__(self, other):
    if isinstance(other, NumMolDescriptorValue):
      return self.value > other.value
    else:
      return self.value > other

  def __lt__(self, other):
    if isinstance(other, NumMolDescriptor):
      return self.value < other.value
    else:
      return self.value < other

class OrdMolDescriptorValue(models.Model):
  '''Contains the ordinal value of a descriptor for a compound'''

  class Meta:
    app_label="DRP"
    verbose_name='Ordinal Molecular Descriptor Value'
    unique_together=('descriptor', 'compound')

  descriptor = models.ForeignKey(BoolMolDescriptor)
  compound = models.ForeignKey(Compound)
  value=models.IntegerField(null=True)

  def __eq__(self, other):
    if isinstance(other, OrdMolDescriptorValue):
      return self.value == other.value 
    else:
      return self.value == other

  def __gt__(self, other):
    if isinstance(other, OrdMolDescriptorValue):
      return self.value > other.value
    else:
      return self.value > other

  def __lt__(self, other):
    if isinstance(other, OrdMolDescriptor):
      return self.value < other.value
    else:
      return self.value < other
