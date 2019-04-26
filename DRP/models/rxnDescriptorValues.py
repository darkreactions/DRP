"""A module containign only the DescriptorValue class."""
from django.db import models
from .descriptorValues import CategoricalDescriptorValue, OrdinalDescriptorValue, BooleanDescriptorValue, NumericDescriptorValue
from .rxnDescriptors import CatRxnDescriptor, NumRxnDescriptor, BoolRxnDescriptor, OrdRxnDescriptor
# Needed to allow for circular dependency.
import DRP.models
import DRP.models.performedReaction
import uuid
from django.contrib.auth.models import User


class RxnDescriptorValueQuerySet(models.query.QuerySet):
    """A queryset which represents a collection of concrete values of a Reaction Descriptor."""

    pass
    # def delete(self):
    # trainingModels = DRP.models.StatsModel.objects.filter(descriptors=self.descriptor, testset__in=dataSets.TestSet.objects.filter(reactions__in=set(v.reaction.performedreaction for v in self)))
    # testModels = DRP.models.StatsModel.objects.filter(descriptors=self.descriptor, trainingset__in=dataSets.TrainingSet.objects.filter(reaction__in=set(v.reaction.performedreaction for v in self)))
    # for model in trainingModels|testModels:
    # model.invalidate()


class RxnDescriptorValueManager(models.Manager):
    """A manager which returns the custom queryset class for Reaction Descriptor Values."""

    def get_queryset(self):
        """Return the correct queryset class."""
        return RxnDescriptorValueQuerySet(self.model, using=self._db)


def rxnUid():
    """Return a unique identifier for a reaction descriptor value."""
    uid = uuid.uuid4()
    while CatRxnDescriptorValue.objects.filter(uid=uid).count() > 0 or BoolRxnDescriptorValue.objects.filter(uid=uid).count() > 0 or NumRxnDescriptorValue.objects.filter(uid=uid).count() > 0 or OrdRxnDescriptorValue.objects.filter(uid=uid).count() > 0:
        uid = uuid.uuid4()
    return uid


class RxnDescriptorValue(models.Model):
    """A class to contain Relationships between Reactions and their descriptors."""

    class Meta:
        app_label = "DRP"
        abstract = True

    uid = models.CharField(max_length=36, default=rxnUid, primary_key=True)
    objects = RxnDescriptorValueManager()
    reaction = models.ForeignKey(
        "DRP.Reaction", unique=False, on_delete=models.PROTECT)

    # def save(self, *args, **kwargs):
    # if self.pk is not None:
    # pass
# try:
#        trainingModels = DRP.models.StatsModel.objects.filter(descriptors=self.descriptor, testset__in=dataSets.TestSet.objects.filter(reactions=self.reaction.performedreaction))
#        testModels = DRP.models.StatsModel.objects.filter(descriptors=self.descriptor, trainingset__in=dataSets.TrainingSet.objects.filter(reaction=self.reaction.performedreaction))
# for model in trainingModels|testModels:
# model.invalidate()
# except pr.PerformedReaction.DoesNotExist:
# pass # fine, we don't care, no need to pass this on.
    # super(RxnDescriptorValue, self).save(*args, **kwargs)

    # def delete(self):
#    trainingModels = DRP.models.StatsModel.objects.filter(descriptors=self.descriptor, testset__in=dataSets.TestSet.objects.filter(reactions=self.reaction.performedreaction))
#    testModels = DRP.models.StatsModel.objects.filter(descriptors=self.descriptor, trainingset__in=dataSets.TrainingSet.objects.filter(reaction=self.reaction.performedreaction))
# for model in trainingModels|testModels:
# model.invalidate()
# model.save()
    # super(RxnDescriptorValue, self).delete()


class CatRxnDescriptorValue(CategoricalDescriptorValue, RxnDescriptorValue):
    """Contains the value of a categorical descriptor for a reaction."""

    descriptorClass = CatRxnDescriptor

    class Meta:
        app_label = "DRP"
        verbose_name = 'Categorical Reaction Descriptor Value'
        unique_together = ('reaction', 'descriptor')


class BoolRxnDescriptorValue(BooleanDescriptorValue, RxnDescriptorValue):
    """Contains the value of a boolean descriptor for a reaction."""

    descriptorClass = BoolRxnDescriptor

    class Meta:
        app_label = "DRP"
        verbose_name = 'Boolean Reaction Descriptor Value'
        unique_together = ('reaction', 'descriptor')


class NumRxnDescriptorValue(NumericDescriptorValue, RxnDescriptorValue):
    """Contains the numeric value of a descriptor for a reaction."""

    descriptorClass = NumRxnDescriptor

    class Meta:
        app_label = "DRP"
        verbose_name = 'Numeric Reaction Descriptor Value'
        unique_together = ('reaction', 'descriptor')


class OrdRxnDescriptorValue(OrdinalDescriptorValue, RxnDescriptorValue):
    """Contains the ordinal value of a descriptor for a reaction."""

    descriptorClass = OrdRxnDescriptor

    class Meta:
        app_label = "DRP"
        verbose_name = 'Ordinal Reaction Descriptor Value'
        unique_together = ('reaction', 'descriptor')

    rater = models.ForeignKey(User, on_delete=models.PROTECT)
