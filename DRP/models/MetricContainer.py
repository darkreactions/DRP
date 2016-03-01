
from django.db import models
from django.conf import settings
from rxnDescriptors import BoolRxnDescriptor, OrdRxnDescriptor, NumRxnDescriptor, CatRxnDescriptor
from dataSets import DataSet
import importlib

metricVisitors = {visitor:importlib.import_module(settings.METRIC_VISITOR_DIR + "." + visitor) for visitor in settings.METRIC_VISITORS}

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

#class TransformedDescriptorAttribute(object):

    ## This style is used to have behavior identical to outcome and regular descriptor properties and the ModelContainer class
    
    #def __get__(self, metricContainer, metricContainerType=None):
        #return metricContainer.transformedRxnDescriptors.all()

    #def __set__(self, metricContainer, descriptors):
        #metricContainer.transformedRxnDescriptors.clear()
        #for descriptor in descriptors:
            #desc = None
            #try:
                #desc = NumRxnDescriptor.objects.get(id=descriptor.id)
                #metricContainer.transformedRxnDescriptors.add(desc)
            #except NumRxnDescriptor.DoesNotExist:
                #pass

            #if desc is None:
                #raise ValueError('An invalid object was assigned as a descriptor')

    #def __delete__(self, metricContainer):
        #metricContainer.transformedRxnDescriptors.clear()


class MetricContainer(models.Model):
    
    class Meta:
        app_label = 'DRP'
        
    description = models.TextField()
    metricVisitor = models.CharField(max_length=255)
    startTime = models.DateTimeField(default=None, null=True)
    endTime = models.DateTimeField(default=None, null=True)
    fileName = models.FileField(upload_to='metrics', max_length=200, blank=True)
    """The filename in which this model is stored"""
    invalid = models.BooleanField(default=False) 
    trainingSet = models.ForeignKey(DataSet, related_name='trainingSetForMetric')
    built = models.BooleanField('Has the build procedure been called with this container?', editable=False, default=False)

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

    #transformedDescriptors = TransformedDescriptorAttribute()
    transformedRxnDescriptors = models.ManyToManyField(NumRxnDescriptor, related_name='transformedByMetric')


    def build(self, predictors, responses, verbose=False):
        if self.built:
            raise RuntimeError("Cannot build a metric that has already been built.")

        self.descriptors = predictors
        self.outcomeDescriptors = responses

        predictorHeaders = [d.csvHeader for d in self.descriptors]
        responseHeaders = [d.csvHeader for d in outcomeDescriptors]

        metricVisitor = metricVisitors[self.metricVisitor].MetricVisitor()
        self.startTime = datetime.datetime.now()
        self.fileName = os.path.join(settings.METRIC_DIR, '{}_{}'.format(self.metricVisitor, self.pk))
        if verbose:
            print "{}, saving to {}, training...".format(self.startTime, self.fileName)
        metricVisitor.train(traningSet.reactions.all(), predictorHeaders, responseHeaders, self.fileName)
        self.endTime = datetime.datetime.now()
        if verbose:
            print "\t...Trained. Finished at {}.".format(self.endTime)
        self.built = True
        

    def transform(self, reactions, verbose=False):
        metricVisitor = metricVisitors[self.metricVisitor].MetricVisitor()
        metricVisitor.recover(self.fileName)
        if verbose:
            print "Transforming..."
        transformed = metricVisitor.transform(reactions)

        if verbose:
            print "\t...transformed"

        if not self.transformedDescriptors.exists():
            if verbose:
                print "No existing descriptors found for this metric container. Creating them..."
            for col in range(transformed.shape[1]):
                desc = NumRxnDescriptor(heading='transform_{}_metric_{}'.format(self.pk, col),
                                        name='{} transformed descriptor for metric container {}'.format(col, self.pk))
                desc.save()
                self.transformedRxnDescriptors.add(desc)
            if verbose:
                print "\t...created"
        elif self.transformedDescriptors.count() != transformed.shape[1]:
            raise RuntimeError("Metric container already has descriptors, but not the same number as the number of columns in data transformed by the metric visitor")
        elif verbose:
            print "Existing transformed descriptors found. Skipping creation."

        if verbose:
            print "Inputting values for given reactions..."
        values = []
        for i, rxn in enumerate(reactions):
            for j, desc in enumerate(self.transformedRxnDescriptors.order_by('pk')):
                val = desc.createValue(rxn, transformed[i,j])
                values.append(val)
        NumRxnDescriptorValue.objects.bulk_create(values)
        if verbose:
            print "\t...finished"
        
