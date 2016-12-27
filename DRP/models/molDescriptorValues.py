"""A module containign only the DescriptorValue class."""
from django.db import models
from .descriptorValues import CategoricalDescriptorValue, BooleanDescriptorValue, NumericDescriptorValue, OrdinalDescriptorValue
# from Compound import DRP.Compound - retain this line for clarity
from django.core.exceptions import ValidationError
import DRP.models
import uuid

#class MolDescriptorValueQuerySet(models.query.QuerySet):
#    """Class to represent a group of Molecular Descriptor Values."""

#    def delete(self, recalculate_reactions=True):
#        """Deletion of a queryset of Descriptor Values should change calculations which arise from them."""
#        if recalculate_reactions:
#            compounds = set(d.compound for d in self)
#        super(MolDescriptorValueQuerySet, self).delete()
#        if recalculate_reactions:
#            for reaction in DRP.models.Reaction.objects.filter(compounds__in=compounds):
#                reaction.save()  # recalculate descriptors
#            for reaction in DRP.models.PerformedReaction.objects.filter(compounds__in=compounds):
#                reaction.save()  # invalidate models


class MolDescriptorValueManager(models.Manager):
    """Manager class to return the custom queryset for MolDescriptorValues."""

    pass
#    def get_queryset(self):
#        """Return the custom queryset."""
#        return MolDescriptorValueQuerySet(self.model, using=self._db)

def molUid():
    """Returns a unique identifier for a molecular descriptor value."""
    uid = uuid.uuid1()
    return uid

class MolDescriptorValue(models.Model):
    """The django model representing the value of a given descriptor for a given compound."""

    class Meta:
        app_label = 'DRP'
        abstract = True


    uid = models.CharField(max_length=36, default=molUid, primary_key=True)
    objects = MolDescriptorValueManager()
    compound = models.ForeignKey('DRP.Compound')

    def delete(self, recalculate_reactions=True):
        """Delete this compound and recalculate each reaction which it is involvedin."""
        reactions = set(self.compound.reaction_set.all())
        super(MolDescriptorValue, self).delete()
        if recalculate_reactions:
            for reaction in reactions:
                reaction.save()  # recalculate descriptors
                try:
                    reaction.performedreaction.save()
                except DRP.models.performedReaction.PerformedReaction.DoesNotExist:
                    pass  # we don't care

    def __str__(self):
        """Return the name, compound and value as the unicode rep."""
        return '{} for {} is {}'.format(self.descriptor.name, str(self.compound), self.value)


class CatMolDescriptorValue(CategoricalDescriptorValue, MolDescriptorValue):
    """The value of a categorical descriptor for a compound."""

    class Meta:
        app_label = "DRP"
        verbose_name = 'Categorical Molecular Descriptor Value'
        unique_together = ('descriptor', 'compound')

    def __str__(self):
        """Return the value of the value for the unicode rep."""
        return self.value.value


class BoolMolDescriptorValue(BooleanDescriptorValue, MolDescriptorValue):
    """The value of a boolean descriptor for a compound."""

    class Meta:
        app_label = "DRP"
        verbose_name = 'Boolean Molecular Descriptor Value'
        unique_together = ('descriptor', 'compound')


class NumMolDescriptorValue(NumericDescriptorValue, MolDescriptorValue):
    """The numeric value of a descriptor for a compound."""

    class Meta:
        app_label = "DRP"
        verbose_name = 'Numeric Molecular Descriptor Value'
        unique_together = ('descriptor', 'compound')


class OrdMolDescriptorValue(OrdinalDescriptorValue, MolDescriptorValue):
    """The ordinal value of a descriptor for a compound."""

    class Meta:
        app_label = "DRP"
        verbose_name = 'Ordinal Molecular Descriptor Value'
        unique_together = ('descriptor', 'compound')
