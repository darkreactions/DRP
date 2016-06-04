'''A module containing only the Reaction class'''
from django.db import models
from LabGroup import LabGroup
from Compound import Compound
from querysets import CsvQuerySet, ArffQuerySet, MultiQuerySet
from descriptors import BooleanDescriptor, NumericDescriptor, CategoricalDescriptor, OrdinalDescriptor
from rxnDescriptorValues import BoolRxnDescriptorValue, NumRxnDescriptorValue, OrdRxnDescriptorValue, CatRxnDescriptorValue
from rxnDescriptors import BoolRxnDescriptor, NumRxnDescriptor, OrdRxnDescriptor, CatRxnDescriptor
from itertools import chain
from CompoundRole import CompoundRole
from collections import OrderedDict
import DRP
import importlib
from django.conf import settings
import gc
from django.db.models.functions import Concat

descriptorPlugins = [importlib.import_module(plugin) for
                     plugin in settings.RXN_DESCRIPTOR_PLUGINS]


class ReactionQuerySet(CsvQuerySet, ArffQuerySet):

    def __init__(self, model=None, **kwargs):
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

    def csvHeaders(self, whitelist=None):
        """Generates the header row information for the CSV"""
        headers = super(ReactionQuerySet, self).csvHeaders(whitelist)
        m = Reaction.objects.all().maxReactantCount()
        for i in range(0, m):
            h = self._getCompoundQuantityHeaderOrder(i)
            if whitelist is None or h in whitelist:
                headers.extend(h)

        return headers

    def arffHeaders(self, whitelist=None):
        """generates headers for the arff file"""
        headers = super(ReactionQuerySet, self).arffHeaders(whitelist)
        m = Reaction.objects.all().maxReactantCount()
        for i in range(0, m):
            compound_label = 'compound_{}'.format(i)
            if whitelist is None or compound_label in whitelist:
                headers[compound_label] = '@attribute {} string'.format(compound_label)
            role_label = 'compound_{}_role'.format(i)
            if whitelist is None or role_label in whitelist:
                headers[role_label] = '@attribute {} {{{}}}'.format(role_label, ','.join(('"{}"'.format(role) for role in CompoundRole.objects.all())))
            amount_label = 'compound_{}_amount'.format(i)
            if whitelist is None or amount_label in whitelist:
                headers[amount_label] = '@attribute {} NUMERIC'.format(amount_label)

        return headers

    def expandedArffHeaders(self, whitelist=None):
        headers = self.arffHeaders(whitelist)
        headers.update(OrderedDict(((d.csvHeader, d.arffHeader) for d in self.descriptors.filter(csvHeader__in=whitelist))))
        return headers

    def expandedCsvHeaders(self, whitelist=None):
        """Generates the expanded header for the csv"""
        if whitelist is not None:
            return self.csvHeaders(whitelist) + [ d.csvHeader for d in self.descriptors.filter(csvHeader__in=whitelist) ]
        else:
            return self.csvHeaders(whitelist) + [ d.csvHeader for d in self.descriptors ]

    #def descriptors(self):
        #"""returns the descriptor which have relationship to the queryset"""
        
        #return chain(
            #BooleanDescriptor.objects.filter(boolrxndescriptorvalue__isnull=False).distinct(),
            #NumericDescriptor.objects.filter(numrxndescriptorvalue__isnull=False).distinct(),
            #OrdinalDescriptor.objects.filter(ordrxndescriptorvalue__isnull=False).distinct(),
            #CategoricalDescriptor.objects.filter(catrxndescriptorvalue__isnull=False).distinct(),
        #)

    @property
    def descriptors(self):
        """returns the descriptor which have relationship to the queryset"""
        return MultiQuerySet(BoolRxnDescriptor.objects.all(),
                                NumRxnDescriptor.objects.all(),
                                OrdRxnDescriptor.objects.all(),
                                CatRxnDescriptor.objects.all()
                                )

    def rows(self, expanded, whitelist=None):
        if expanded:
            reactions = self
            if whitelist is not None:
                # This whole nonsense just grabs the descriptor values we actually want. Let's break it down
                # First we annotate each descriptor value with its descriptor's csvHeader using the Concat operation
                # Then we filter on the csvHeader so it's in the whitelist
                # Then we prefetch descriptor values matching this queryset for each reaction and assign this to filtered __vals
                # Then we can just iterate through these values to get what we actually want!
                qs = BoolRxnDescriptorValue.objects.annotate(descCsvHeader=Concat('descriptor__heading', models.Value('_'), 'descriptor__calculatorSoftware', models.Value('_'), 'descriptor__calculatorSoftwareVersion')).filter(descCsvHeader__in=whitelist)
                reactions = reactions.prefetch_related(models.Prefetch('boolrxndescriptorvalue_set', queryset=qs, to_attr='filtered_boolvals'))
                
                qs = NumRxnDescriptorValue.objects.annotate(descCsvHeader=Concat('descriptor__heading', models.Value('_'), 'descriptor__calculatorSoftware', models.Value('_'), 'descriptor__calculatorSoftwareVersion')).filter(descCsvHeader__in=whitelist)
                reactions = reactions.prefetch_related(models.Prefetch('numrxndescriptorvalue_set', queryset=qs, to_attr='filtered_numvals'))
                
                qs = OrdRxnDescriptorValue.objects.annotate(descCsvHeader=Concat('descriptor__heading', models.Value('_'), 'descriptor__calculatorSoftware', models.Value('_'), 'descriptor__calculatorSoftwareVersion')).filter(descCsvHeader__in=whitelist)
                reactions = reactions.prefetch_related(models.Prefetch('ordrxndescriptorvalue_set', queryset=qs, to_attr='filtered_ordvals'))
                
                qs = CatRxnDescriptorValue.objects.annotate(descCsvHeader=Concat('descriptor__heading', models.Value('_'), 'descriptor__calculatorSoftware', models.Value('_'), 'descriptor__calculatorSoftwareVersion')).filter(descCsvHeader__in=whitelist)
                reactions = reactions.prefetch_related(models.Prefetch('catrxndescriptorvalue_set', queryset=qs, to_attr='filtered_catvals'))
            else: 
                reactions = reactions.prefetch_related('boolrxndescriptorvalue_set__descriptor')
                reactions = reactions.prefetch_related('catrxndescriptorvalue_set__descriptor')
                reactions = reactions.prefetch_related('ordrxndescriptorvalue_set__descriptor')
                reactions = reactions.prefetch_related('numrxndescriptorvalue_set__descriptor')
            reactions = reactions.prefetch_related('compounds')
            
            for item in reactions.batch_iterator():
                row = {field.name:getattr(item, field.name) for field in self.model._meta.fields}
                if whitelist is not None:
                    row.update({dv.descCsvHeader:dv.value for dv in item.filtered_boolvals})
                    row.update({dv.descCsvHeader:dv.value for dv in item.filtered_numvals})
                    row.update({dv.descCsvHeader:dv.value for dv in item.filtered_ordvals})
                    row.update({dv.descCsvHeader:dv.value for dv in item.filtered_catvals})
                    i=0
                    for compound in item.compounds.all():
                        compound_num = 'compound_{}'.format(i)
                        if compound_num in whitelist:
                            row[compound_num] = compound.name
                        i+=1
                    yield row
                else:
                    row.update({dv.descriptor.csvHeader:dv.value for dv in item.descriptorValues})
                    i=0
                    for compound in item.compounds.all():
                        row['compound_{}'.format(i)] = compound.name
                        i+=1
                    yield row
        else:
            for row in super(ReactionQuerySet, self).rows(expanded):
                yield row



    # From https://djangosnippets.org/snippets/1949/
    def batch_iterator(self, chunksize=5000):
        '''
        Iterate over a Django Queryset ordered by the primary key
    
        This method loads a maximum of chunksize (default: 5000) rows in it's
        memory at the same time while django normally would load all rows in it's
        memory. Using the iterator() method only causes it to not preload all the
        classes.
    
        Note that the implementation of the iterator does not support ordered query sets.
        '''
        pk = 0
        last_pk = self.order_by('-pk')[0].pk
        queryset = self.order_by('pk')
        while pk < last_pk:
            for row in queryset.filter(pk__gt=pk)[:chunksize]:
                pk = row.pk
                yield row
            gc.collect()

    def calculate_descriptors(self, verbose=False, plugins=None, **kwargs):
        if verbose:
            print "Calculating descriptors for {} reactions".format(self.count())
        for plugin in descriptorPlugins:
            if plugins is None or plugin.__name__ in plugins:
                if verbose:
                    print "Calculating for plugin: {}".format(plugin)
                plugin.calculate_many(self, verbose=verbose, **kwargs)
                if verbose:
                    print "Done with plugin: {}\n".format(plugin)

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
    compounds = models.ManyToManyField(Compound, through="CompoundQuantity")
    notes = models.TextField(blank=True)
    labGroup = models.ForeignKey(LabGroup)
    calcDescriptors = True #this is to cope with a hideous problem in xml serialization in the management commands


    def save(self, calcDescriptors=False, *args, **kwargs):
        super(Reaction, self).save(*args, **kwargs)
        if calcDescriptors and self.calcDescriptors:
            for plugin in descriptorPlugins:
                plugin.calculate(self)
                
    #@property
    #def descriptorValues(self):
        #return chain(self.boolrxndescriptorvalue_set.all(), self.numrxndescriptorvalue_set.all(), self.ordrxndescriptorvalue_set.all(), self.catrxndescriptorvalue_set.all())
        
    @property
    def descriptorValues(self):
        return MultiQuerySet(self.boolrxndescriptorvalue_set.all(), self.numrxndescriptorvalue_set.all(), self.ordrxndescriptorvalue_set.all(), self.catrxndescriptorvalue_set.all())

    def __unicode__(self):
        return "Reaction_{}".format(self.id)
