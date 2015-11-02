'''A module containing only the Reaction class'''
from django.db import models
from LabGroup import LabGroup
from Compound import Compound
from querysets import CsvQuerySet, ArffQuerySet
from descriptors import BooleanDescriptor, NumericDescriptor, CategoricalDescriptor, OrdinalDescriptor
from rxnDescriptorValues import BoolRxnDescriptorValue, NumRxnDescriptorValue, OrdRxnDescriptorValue, CatRxnDescriptorValue
from itertools import chain
from CompoundRole import CompoundRole
from collections import OrderedDict

import importlib
from django.conf import settings

descriptorPlugins = [importlib.import_module(plugin) for
                     plugin in settings.RXN_DESCRIPTOR_PLUGINS]


class ReactionQuerySet(CsvQuerySet, ArffQuerySet):

    def __init__(self, model = None, **kwargs):
        """Initialises the queryset"""
        model = Reaction if model is None else model
        super(ReactionQuerySet, self).__init__(model=model, **kwargs)

    def maxReactantCount(self):
        """Gives a count of the maximum number of reactions associated with this queryset"""
        m = self.annotate(compoundQuantityCount=models.Count('compoundquantity')).aggregate(max=models.Max('compoundQuantityCount'))['max']
        if m is None:
            return 0
        return m

    def _getCompoundQuantityHeaderOrder(self, i):
        return ['compound_{}'.format(i), 'compound_{}_role'.format(i), 'compound_{}_amount'.format(i)]

    @property
    def csvHeaders(self):
        """Generates the header row information for the CSV"""
        headers = super(ReactionQuerySet, self).csvHeaders
        m = Reaction.objects.all().maxReactantCount()
        for i in range(0, m):
            headers.extend(self._getCompoundQuantityHeaderOrder(i))

        return headers

    @property
    def arffHeaders(self):
        """generates headers for the arff file"""
        headers = super(ReactionQuerySet, self).arffHeaders
        m = Reaction.objects.all().maxReactantCount()
        for i in range(0, m):
            compound_label = 'compound_{}'.format(i)
            headers[compound_label] = '@attribute {} string'.format(compound_label)
            role_label = 'compound_{}_role'.format(i)
            headers[role_label] = '@attribute {} {{{}}}'.format(role_label, ','.join(('"{}"'.format(role) for role in CompoundRole.objects.all())))
            amount_label = 'compound_{}_amount'.format(i)
            headers[amount_label] = '@attribute {} NUMERIC'.format(amount_label)

        return headers

    @property
    def expandedArffHeaders(self):
        headers = self.arffHeaders
        headers.update(OrderedDict(((d.csvHeader, d.arffHeader) for d in self.descriptors())))
        return headers

    @property
    def expandedCsvHeaders(self):
        """Generates the expanded header for the csv"""
        return self.csvHeaders + [ d.csvHeader for d in self.descriptors() ]

    def descriptors(self):
        """returns the descriptor which have relationship to the queryset"""
        return chain(BooleanDescriptor.objects.filter(boolrxndescriptorvalue__in=BoolRxnDescriptorValue.objects.filter(reaction__in=self)), NumericDescriptor.objects.filter(numrxndescriptorvalue__in=NumRxnDescriptorValue.objects.filter(reaction__in=self)), OrdinalDescriptor.objects.filter(ordrxndescriptorvalue__in=OrdRxnDescriptorValue.objects.filter(reaction__in=self)), CategoricalDescriptor.objects.filter(catrxndescriptorvalue__in=CatRxnDescriptorValue.objects.filter(reaction__in=self)))


class ReactionManager(models.Manager):
    """A custom manager for the Reaction Class which permits the creation of entries to and from CSVs"""
    use_for_related_fields = True

    def get_queryset(self):
        return ReactionQuerySet()


class Reaction(models.Model):
  '''A base class on which PerformedReactions and RecommendedReactions are built,
  contains common information to each in a table with an automatically
  generated one to one relationship with the subclasses.
  '''

  class Meta:
    app_label="DRP"

  objects = ReactionManager()
  compounds=models.ManyToManyField(Compound, through="CompoundQuantity")
  notes=models.TextField(blank=True)
  labGroup=models.ForeignKey(LabGroup)

  def save(self, *args, **kwargs):
    super(Reaction, self).save(*args, **kwargs)
    for plugin in descriptorPlugins:
      plugin.calculate(self)
