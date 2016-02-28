
from django.db import models
from rxnDescriptors import BoolRxnDescriptor, OrdRxnDescriptor, NumRxnDescriptor, CatRxnDescriptor


class DescriptorAttribute(object):

    def __get__(self, metricContainer):
        return chain(metricContainer.boolRxnDescriptors.all(), metricContainer.ordRxnDescriptors.all(), metricContainer.catRxnDescriptors.all(), metricContainer.numRxnDescriptors.all())

    def __set__(self, metricContainer, descriptors):
        metricContainer.boolRxnDescriptors.clear()
        metricContainer.ordRxnDescriptors.clear()
        metricContainer.catRxnDescriptors.clear()
        metricContainer.numRxnDescriptors.clear()
        for descriptor in descriptors:
            desc = None
            try:
                desc = BoolRxnDescriptor.objects.get(id=descriptor.id)
                metricContainer.boolRxnDescriptors.add(desc)
            except BoolRxnDescriptor.DoesNotExist:
                pass

            try:
                desc = OrdRxnDescriptor.objects.get(id=descriptor.id)
                metricContainer.ordRxnDescriptors.add(desc)
            except OrdRxnDescriptor.DoesNotExist:
                pass

            try:
                desc = CatRxnDescriptor.objects.get(id=descriptor.id)
                metricContainer.catRxnDescriptors.add(desc)
            except CatRxnDescriptor.DoesNotExist:
                pass

            try:
                desc = NumRxnDescriptor.objects.get(id=descriptor.id)
                metricContainer.numRxnDescriptors.add(desc)
            except NumRxnDescriptor.DoesNotExist:
                pass

            if desc is None:
                print descriptor.heading
                print type(descriptor)
                raise ValueError('An invalid object was assigned as a descriptor')

    def __delete__(self, metricContainer):
        metricContainer.boolRxnDescriptors.clear()
        metricContainer.numRxnDescriptors.clear()
        metricContainer.catRxnDescriptors.clear()
        metricContainer.ordRxnDescriptors.clear()

class OutcomeDescriptorAttribute(object):

    def __get__(self, metricContainer, metricContainerType=None):
        return chain(metricContainer.outcomeBoolRxnDescriptors.all(), metricContainer.outcomeOrdRxnDescriptors.all(), metricContainer.outcomeCatRxnDescriptors.all(), metricContainer.outcomeNumRxnDescriptors.all())

    def __set__(self, metricContainer, descriptors):
        metricContainer.outcomeBoolRxnDescriptors.clear()
        metricContainer.outcomeOrdRxnDescriptors.clear()
        metricContainer.outcomeCatRxnDescriptors.clear()
        metricContainer.outcomeNumRxnDescriptors.clear()
        for descriptor in descriptors:
            desc = None
            try:
                desc = BoolRxnDescriptor.objects.get(id=descriptor.id)
                metricContainer.outcomeBoolRxnDescriptors.add(desc)
            except BoolRxnDescriptor.DoesNotExist:
                pass

            try:
                desc = OrdRxnDescriptor.objects.get(id=descriptor.id)
                metricContainer.outcomeOrdRxnDescriptors.add(desc)
            except OrdRxnDescriptor.DoesNotExist:
                pass

            try:
                desc = CatRxnDescriptor.objects.get(id=descriptor.id)
                metricContainer.outcomeCatRxnDescriptors.add(desc)
            except CatRxnDescriptor.DoesNotExist:
                pass

            try:
                desc = NumRxnDescriptor.objects.get(id=descriptor.id)
                metricContainer.outcomeNumRxnDescriptors.add(desc)
            except NumRxnDescriptor.DoesNotExist:
                pass

            if desc is None:
                raise ValueError('An invalid object was assigned as a descriptor')

    def __delete__(self, metricContainer):
        metricContainer.outcomeBoolRxnDescriptors.clear()
        metricContainer.outcomeNumRxnDescriptors.clear()
        metricContainer.outcomeCatRxnDescriptors.clear()
        metricContainer.outcomeOrdRxnDescriptors.clear()

class TransformedDescriptorAttribute(object):

    def __get__(self, metricContainer, metricContainerType=None):
        return chain(metricContainer.transformedRxnDescriptors.all())

    def __set__(self, metricContainer, descriptors):
        metricContainer.transformedRxnDescriptors.clear()
        for descriptor in descriptors:
            desc = None
            try:
                desc = NumRxnDescriptor.objects.get(id=descriptor.id)
                metricContainer.transformedRxnDescriptors.add(desc)
            except NumRxnDescriptor.DoesNotExist:
                pass

            if desc is None:
                raise ValueError('An invalid object was assigned as a descriptor')

    def __delete__(self, metricContainer):
        metricContainer.transformedRxnDescriptors.clear()


class MetricContainer(models.Model):
    
    class Meta:
        app_label = 'DRP'
        
    description = models.TextField()
    metricVisitor = models.CharField(max_length=255)
    splitter = models.CharField(max_length=255)
    startTime = models.DateTimeField(default=None, null=True)
    endTime = models.DateTimeField(default=None, null=True)

    trainingSets = models.ForeignKey(DataSet, related_name='trainingSetFor')
    testSets = models.ManyToManyField(DataSet, related_name='testSetsFor')

    descriptors = DescriptorAttribute()
    boolRxnDescriptors = models.ManyToManyField(BoolRxnDescriptor)
    ordRxnDescriptors = models.ManyToManyField(OrdRxnDescriptor)
    catRxnDescriptors = models.ManyToManyField(CatRxnDescriptor)
    numRxnDescriptors = models.ManyToManyField(NumRxnDescriptor)

    outcomeDescriptors = OutcomeDescriptorAttribute()
    outcomeBoolRxnDescriptors = models.ManyToManyField(BoolRxnDescriptor, related_name='outcomeForMetrics')
    outcomeOrdRxnDescriptors = models.ManyToManyField(OrdRxnDescriptor, related_name='outcomeForMetrics')
    outcomeCatRxnDescriptors = models.ManyToManyField(CatRxnDescriptor, related_name='outcomeForMetrics')
    outcomeNumRxnDescriptors = models.ManyToManyField(NumRxnDescriptor, related_name='outcomeForMetrics')

    transformedDescriptors = TransformedDescriptorAttribute()
    transformedRxnDescriptors = models.ManyToManyField(NumRxnDescriptor, related_name='transformedByMetric')
    
    @classmethod
    def create(cls, metricVisitor, description=None, splitter=None, reactions=None, trainingSet=None, testSet=None):
        pass


    def clean(self):
        if (self.splitter is None) ^ (self.reactions is None): # if these are not the same, there's a problem
            raise ValidationError('A full set of reactions must be supplied with a splitter', 'argument_mismatch')
        if not ((self.splitter is None) ^ (self.trainingSets is None)): # if these are not different, there's a problem
            raise ValidationError('Either a splitter or a training set should be provided.', 'argument_mismatch')

    
