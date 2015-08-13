'''A module containign only the DescriptorValue class'''
from django.db import models
from MolDescriptor import CatMolDescriptor, BoolMolDescriptor, NumMolDescriptor, OrdMolDescriptor,CatMolDescriptorPermitted
from Reaction import Reaction
from Compound import Compound
from StatsModel import StatsModel

class CatMolDescriptorValue(models.Model):
  '''Contains the value of a categorical descriptor for a compound'''

  class Meta:
    app_label="DRP"
    verbose_name='Categorical Molecular Descriptor Value'
    unique_together=('descriptor', 'compound')

  descriptor = models.ForeignKey(CatMolDescriptor)
  compound = models.ForeignKey(Compound)
  value = models.ForeignKey(CatMolDescriptorPermitted) #NOTE that the correct permission cannot be enforced at the model or database level! make sure this is done in all forms!

class BoolMolDescriptorValue(models.Model):
  '''Contains the value of a boolean descriptor for a compound'''

  class Meta:
    app_label="DRP"
    verbose_name='Boolean Molecular Descriptor Value'
    unique_together=('descriptor', 'compound')

  descriptor = models.ForeignKey(BoolMolDescriptor)
  compound = models.ForeignKey(Compound)
  value= models.NullBooleanField('Value if descriptor is a boolean', null=True)

class NumMolDescriptorValue(models.Model):
  '''Contains the numeric value of a descriptor for a compound'''

  class Meta:
    app_label="DRP"
    verbose_name='Boolean Molecular Descriptor Value'
    unique_together=('descriptor', 'compound')

  descriptor = models.ForeignKey(BoolMolDescriptor)
  compound = models.ForeignKey(Compound)
  value=models.FloatField(null=True)

class OrdMolDescriptorValue(models.Model):
  '''Contains the ordinal value of a descriptor for a compound'''

  class Meta:
    app_label="DRP"
    verbose_name='Ordinal Molecular Descriptor Value'
    unique_together=('descriptor', 'compound')

  descriptor = models.ForeignKey(BoolMolDescriptor)
  compound = models.ForeignKey(Compound)
  value=models.IntegerField(null=True)

