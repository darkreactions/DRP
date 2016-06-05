from django.db import models


class RxnDescriptorSet(models.model):

    """A class for describing a group of descriptors."""

    class Meta:
        app_label = 'DRP'

    description = models.TextField(blank=True, null=False)
    boolRxnDescriptors = models.ManyToManyField(BoolRxnDescriptor)
    ordRxnDescriptors = models.ManyToManyField(OrdRxnDescriptor)
    catRxnDescriptors = models.ManyToManyField(CatRxnDescriptor)
    numRxnDescriptors = models.ManyToManyField(NumRxnDescriptor)

    def add(self, desc):
        desc = None
        try:
            desc = BoolRxnDescriptor.objects.get(id=descriptor.id)
            self.boolRxnDescriptors.add(desc)
        except BoolRxnDescriptor.DoesNotExist:
            pass

        try:
            desc = OrdRxnDescriptor.objects.get(id=descriptor.id)
            self.ordRxnDescriptors.add(desc)
        except OrdRxnDescriptor.DoesNotExist:
            pass

        try:
            desc = CatRxnDescriptor.objects.get(id=descriptor.id)
            self.catRxnDescriptors.add(desc)
        except CatRxnDescriptor.DoesNotExist:
            pass

        try:
            desc = NumRxnDescriptor.objects.get(id=descriptor.id)
            self.numRxnDescriptors.add(desc)
        except NumRxnDescriptor.DoesNotExist:
            pass

        if desc is None:
            raise ValueError('An invalid object was assigned as a descriptor')

    def set(self, descriptors):
        self.boolRxnDescriptors.clear()
        self.ordRxnDescriptors.clear()
        self.catRxnDescriptors.clear()
        self.numRxnDescriptors.clear()
        for descriptor in descriptors:
            self.add(descriptor)

    def get(self):
        return chain(self.boolRxnDescriptors.all(), self.ordRxnDescriptors.all(), self.catRxnDescriptors.all(), self.numRxnDescriptors.all())
