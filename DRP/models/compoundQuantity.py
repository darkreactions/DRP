"""Module containing only the CompoundQuantities Class."""
from django.db import models
from django import forms
from django.core.exceptions import ValidationError
from .compound import Compound
from .compoundRole import CompoundRole
from .reaction import Reaction
from .performedReaction import PerformedReaction
from .validators import GreaterThanValidator


class CompoundQuantityQuerySet(models.query.QuerySet):

    """A queryset representing a collection of compounds quantities."""

    def delete(self):
        """Force the re-save of reactions pertinent to these compound quantities on deletion."""
        reactions = Reaction.objects.filter(compoundquantity_set__in=self)
        for reaction in reactions:
            reaction.save()  # recalculate descriptors
            try:
                reaction.performedreaction.save()
            except PerformedReaction.DoesNotExist:
                pass  # we don't care about this outcome
        super(CompdoundQuantityQuerySet, self).delete()


class CompoundQuantityManager(models.Manager):

    """A manager for CompoundQuantity objects."""

    def get_queryset(self):
        """Return the appropriate custom queryset."""
        return CompoundQuantityQuerySet(self.model, using=self._db)


class CompoundQuantity(models.Model):

    """
    A class to contain the relationship between a reaction and a compound.

    Contains the amount of a given compound used in a reaction
    with the applicable units. At present, no unit convention is enforced.
    """

    class Meta:
        app_label = 'DRP'
        unique_together = ('reaction', 'role', 'amount')

    compound = models.ForeignKey(Compound, on_delete=models.PROTECT)
    reaction = models.ForeignKey(Reaction)
    role = models.ForeignKey(CompoundRole)
    amount = models.DecimalField(null=True, blank=True, max_digits=12, decimal_places=5,
                                 help_text="(in mmoles, 5 decimal places)", validators=[GreaterThanValidator(0)])

    def save(self, calcDescriptors=False, invalidate_models=True, *args, **kwargs):
        """Re-save associated reactions dependent upon this quantity as this will cause descriptor values to change."""
        super(CompoundQuantity, self).save(*args, **kwargs)
        try:
            self.reaction.performedreaction.save(
                calcDescriptors=calcDescriptors, invalidate_models=invalidate_models)  # invalidate models
        except PerformedReaction.DoesNotExist:
            # descriptor recalculation
            self.reaction.save(calcDescriptors=calcDescriptors)

    def delete(self):
        """Re-save associated reactions dependent upon this quantity as this will cause descriptor values to change."""
        try:
            self.reaction.performedreaction.save()  # invalidate models
        except PerformedReaction.DoesNotExist:
            self.reaction.save()  # descriptor recalculation
        super(CompoundQuantity, self).save()

    def __unicode__(self):
        """Return the compound, amount and reaction as a unicode representation."""
        return u'{} {} with role {} in {}'.format(self.amount, self.compound, self.role, self.reaction)
