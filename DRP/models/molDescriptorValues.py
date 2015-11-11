'''A module containign only the DescriptorValue class'''
from django.db import models
from descriptorValues import CategoricalDescriptorValue, BooleanDescriptorValue, NumericDescriptorValue, OrdinalDescriptorValue
#from Compound import DRP.Compound - retain this line for clarity
from django.core.exceptions import ValidationError
import StatsModel
import PerformedReaction
import DRP.models

class MolDescriptorValueQuerySet(models.query.QuerySet):

  def delete(self):
    compounds = set(d.compound for d in self)
    for reaction in DRP.models.Reaction.objects.filter(compounds__in=compounds):
      reaction.save() #recalculate descriptors
    for reaction in DRP.models.PerformedReaction.objects.filter(compounds__in=compounds):
      reaction.save() #invalidate models

class MolDescriptorValueManager(models.Manager):

  def get_queryset(self):
    return MolDescriptorValueQuerySet(self.model, using=self._db)

class MolDescriptorValue(models.Model):

  class Meta:
    app_label ='DRP'
    abstract=True

  objects = MolDescriptorValueManager()
  compound = models.ForeignKey('DRP.Compound')

  def delete(self):
    for reaction in self.compound.reaction_set.all():
      reaction.save() #recalculate descriptors
      try:
        reaction.performedreaction.save()
      except PerformedReaction.PerformedReaction.DoesNotExist:
        pass #we don't care

  def __str__(self):
    return '{} for {} is {}'.format(self.descriptor.name, str(self.compound), self.value)

class CatMolDescriptorValue(CategoricalDescriptorValue, MolDescriptorValue):
  '''Contains the value of a categorical descriptor for a compound'''

  class Meta:
    app_label="DRP"
    verbose_name='Categorical Molecular Descriptor Value'
    unique_together=('descriptor', 'compound')

  def __unicode__(self):
    return self.value.value

class BoolMolDescriptorValue(BooleanDescriptorValue, MolDescriptorValue):
  '''Contains the value of a boolean descriptor for a compound'''

  class Meta:
    app_label="DRP"
    verbose_name='Boolean Molecular Descriptor Value'
    unique_together=('descriptor', 'compound')

class NumMolDescriptorValue(NumericDescriptorValue, MolDescriptorValue):
  '''Contains the numeric value of a descriptor for a compound'''

  class Meta:
    app_label="DRP"
    verbose_name='Numeric Molecular Descriptor Value'
    unique_together=('descriptor', 'compound')

class OrdMolDescriptorValue(OrdinalDescriptorValue, MolDescriptorValue):
  '''Contains the ordinal value of a descriptor for a compound'''

  class Meta:
    app_label="DRP"
    verbose_name='Ordinal Molecular Descriptor Value'
    unique_together=('descriptor', 'compound')

