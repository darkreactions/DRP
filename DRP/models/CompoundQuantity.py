'''Module containing only the CompoundQuantities Class'''
from django.db import models
from django import forms
from django.core.exceptions import ValidationError
from Compound import Compound
from CompoundRole import CompoundRole
from Reaction import Reaction
from PerformedReaction import PerformedReaction

class CompoundQuantityQuerySet(models.query.QuerySet):

    def delete(self):
        reactions = Reaction.objects.filter(compoundquantity_set__in=self)
        for reaction in reactions:
            reaction.save() #recalculate descriptors
            try:
                reaction.performedreaction.save()
            except PerformedReaction.DoesNotExist:
                pass #we don't care about this outcome

class CompoundQuantityManager(models.Manager):

    def get_queryset(self):
        return CompoundQuantityQuerySet(self.model, using=self._db)

class CompoundQuantity(models.Model):
    '''A class to contain the relationship between a reaction and a compound,
    and thus to contain the amount of a given compound used in a reaction
    with the applicable units. At present, no unit convention is enforced.
    '''

    class Meta:
        app_label='DRP'
        unique_together=('reaction', 'role', 'amount')

    compound=models.ForeignKey(Compound, on_delete=models.PROTECT)
    reaction=models.ForeignKey(Reaction)
    role=models.ForeignKey(CompoundRole)
    amount=models.FloatField(null=True, blank=True)

    def save(self, *args, **kwargs):
        self.reaction.save() #descriptor recalculation
        try:
            self.reaction.performedreaction.save() #invalidate models
        except PerformedReaction.DoesNotExist:
            pass #we don't care that it doesn't exist
        super(CompoundQuantity, self).save()

    def delete(self):
        self.reaction.save() #descriptor recalculation
        try:
            self.reaction.performedreaction.save() #invalidate models
        except PerformedReaction.DoesNotExist:
            pass #we don't care
        super(CompoundQuantity, self).save()
