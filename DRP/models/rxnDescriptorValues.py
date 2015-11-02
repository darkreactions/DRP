'''A module containign only the DescriptorValue class'''
from django.db import models
from Reaction import Reaction
from PerformedReaction import PerformedReaction
from descriptorValues import CategoricalDescriptorValue, OrdinalDescriptorValue,BooleanDescriptorValue, NumericDescriptorValue
from StatsModel import StatsModel
from itertools import chain

class RxnDescriptorValueQuerySet(models.query.QuerySet):

  def delete(self):
    for model in StatsModel.objects.filter(trainingset__in=self) | StatsModel.objects.filter(testset_in=self):
      model.invalidate()

class RxnDescriptorValueManager(models.Manager):

  def get_queryset(self):
    return RxnDescriptorValueQuerySet(self.model, using=self._db)

class RxnDescriptorValue(models.Model):
  '''Contains Relationships between Reactions and their descriptors'''

  class Meta:
    app_label="DRP"
    abstract=True

  objects = RxnDescriptorValueManager()
  reaction = models.ForeignKey(Reaction, unique=False)
  model=models.ForeignKey(StatsModel, unique=False, null=True, default=None)
  '''If this value was predicted by a statistical model, reference that model'''

  def save(self, *args, **kwargs):
    if self.pk is not None:
      try:
        for model in chain(self.reaction.performedreaction.inTrainingSetFor.all(), self.reaction.performedreaction.inTestSetFor.all()):
          model.invalidate()
      except PerformedReaction.DoesNotExist:
        pass # fine, we don't care, no need to pass this on.
    super(RxnDescriptorValue, self).save(*args, **kwargs)

  def delete(self):
    test = StatsModel.objects.filter(testset__reactions__in=[self])
    train = StatsModel.objects.filter(trainingset__reaction=self)
    for model in chain(test, train):
      model.invalidate()
      model.save()
    super(RxnDescriptorValue, self).delete()

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
