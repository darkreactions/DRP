from descriptors import CategoricalDescriptorPermittedValue, CategoricalDescriptor, NumericDescriptor, BooleanDescriptor, OrdinalDescriptor
from django.core.exceptions import ValidationError
from django.db import models

class CategoricalDescriptorValue(models.Model):

  descriptor = models.ForeignKey(CategoricalDescriptor)
  value = models.ForeignKey(CategoricalDescriptorPermittedValue, null=True)

  class Meta:
    app_label='DRP'
    abstract=True

  def __eq__(self, other):
    if isinstance(other, CategoricalDescriptorValue):
      return self.value.value == other.value.value
    elif isinstance(other, CategoricalDescriptorPermittedValue):
      return self.value.value == other.value
    else:
      return self.value.value == other

  def clean(self):
    if self.value is not None and self.value not in self.descriptor.permittedValues.all():
      raise ValidationError('Invalid Category Described for this Categorical Descriptor', 'invalid_category')

  def save(self, *args, **kwargs):
    self.clean()
    super(CategoricalDescriptorValue, self).save(*args, **kwargs)

class BooleanDescriptorValue(models.Model):

  class Meta:
    app_label='DRP'
    abstract=True

  value= models.NullBooleanField('Value for descriptor', null=True)
  descriptor = models.ForeignKey(BooleanDescriptor)

  def __nonzero__(self):
    return self.value

class NumericDescriptorValue(models.Model):
  
  class Meta:
    app_label='DRP'
    abstract=True

  value=models.FloatField(null=True)
  descriptor=models.ForeignKey(NumericDescriptor)

  def clean(self):
    if self.value is not None:
      if self.descriptor.maximum is not None and self > self.descriptor.maximum:
        raise ValidationError('The provided value is higher than the descriptor maximum', 'value_too_high')
      if self.descriptor.minimum is not None and self < self.descriptor.minimum:
        raise ValidationError('The provided value is lower than the descriptor minimum', 'value_too_low')

  def save(self, *args, **kwargs):
    self.clean()
    super(NumericDescriptorValue, self).save(*args, **kwargs)

  def __nonzero__(self):
    return bool(self.value)

  def __eq__(self, other):
    if isinstance(other, NumericDescriptorValue):
      return self.value == other.value 
    else:
      return self.value == other

  def __gt__(self, other):
    if isinstance(other, NumericDescriptorValue):
      return self.value > other.value
    else:
      return self.value > other

  def __lt__(self, other):
    if isinstance(other, NumericDescriptorValue):
      return self.value < other.value
    else:
      return self.value < other

class OrdinalDescriptorValue(models.Model):

  class Meta:
    app_label='DRP'
    abstract=True

  value=models.IntegerField(null=True)
  descriptor = models.ForeignKey(OrdinalDescriptor)
  
  def clean(self):
    if self.value is not None:
      if self.descriptor.maximum is not None and self > self.descriptor.maximum:
        raise ValidationError('The provided value is higher than the descriptor maximum', 'value_too_high')
      if self.descriptor.minimum is not None and self < self.descriptor.minimum:
        raise ValidationError('The provided value is lower than the descriptor minimum', 'value_too_low')

  def save(self, *args, **kwargs):
    self.clean()
    super(OrdinalDescriptorValue, self).save(*args, **kwargs)

  def __eq__(self, other):
    if isinstance(other, OrdinalDescriptorValue):
      return self.value == other.value 
    else:
      return self.value == other

  def __gt__(self, other):
    if isinstance(other, OrdinalDescriptorValue):
      return self.value > other.value
    else:
      return self.value > other

  def __lt__(self, other):
    if isinstance(other, OrdinalDescriptorValue):
      return self.value < other.value
    else:
      return self.value < other
