'''A module containign only the DescriptorValue class'''
from django.db import models
from MolDescriptor import CatMolDescriptor, BoolMolDescriptor, NumMolDescriptor, OrdMolDescriptor,CatMolDescriptorPermitted
from Reaction import Reaction
from Compound import Compound
from StatsModel import StatsModel
from django.core.exceptions import ValidationError

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
  value = models.ForeignKey(CatMolDescriptorPermitted, null=True)

  def __eq__(self, other):
    if isinstance(other, CatMolDescriptorValue):
      return self.value.value == other.value.value
    elif isinstance(other, CatMolDescriptorPermitted):
      return self.value.value == other.value
    else:
      return self.value.value == other

  def clean(self):
    if self.value not in self.descriptor.permittedValues.all() and self.value is not None:
      raise ValidationError('Invalid Category Described for this Categorical Descriptor', 'invalid_category')

  def save(self, *args, **kwargs):
    self.clean()
    super(CatMolDescriptorValue, self).save(*args, **kwargs)


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

  descriptor = models.ForeignKey(NumMolDescriptor)
  compound = models.ForeignKey(Compound)
  value=models.FloatField(null=True)

  def clean(self):
    if self.value is not None:
      if self.descriptor.maximum is not None and self > self.descriptor.maximum:
        raise ValidationError('The provided value is higher than the descriptor maximum', 'value_too_high')
      if self.descriptor.minimum is not None and self < self.descriptor.minimum:
        raise ValidationError('The provided value is lower than the descriptor minimum', 'value_too_low')

  def save(self, *args, **kwargs):
    self.clean()
    super(NumMolDescriptorValue, self).save(*args, **kwargs)

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

  descriptor = models.ForeignKey(OrdMolDescriptor)
  compound = models.ForeignKey(Compound)
  value=models.IntegerField(null=True)
  
  def clean(self):
    if self.value is not None:
      if self.descriptor.maximum is not None and self > self.descriptor.maximum:
        raise ValidationError('The provided value is higher than the descriptor maximum', 'value_too_high')
      if self.descriptor.minimum is not None and self < self.descriptor.minimum:
        raise ValidationError('The provided value is lower than the descriptor minimum', 'value_too_low')

  def save(self, *args, **kwargs):
    self.clean()
    super(OrdMolDescriptorValue, self).save(*args, **kwargs)

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
