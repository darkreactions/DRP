import django
django.setup()
from DRP.models import (DataSet, DataSetRelation, ModelContainer, MetricContainer, FeatureSelectionContainer, StatsModel,
                        CatRxnDescriptorValue, OrdRxnDescriptorValue, BoolRxnDescriptorValue, NumRxnDescriptorValue,
                        CatMolDescriptorValue, OrdMolDescriptorValue, BoolMolDescriptorValue, NumMolDescriptorValue,
                        CatRxnDescriptor, OrdRxnDescriptor, BoolRxnDescriptor, NumRxnDescriptor,
                        CatMolDescriptor, OrdMolDescriptor, BoolMolDescriptor, NumMolDescriptor,
                        RxnDescriptorValue,
                        Descriptor)


if __name__ == '__main__':
    print "Deleting DataSetRelation"
    DataSetRelation.objects.all().delete()
    print "Deleting DataSet"
    DataSet.objects.all().delete()
    print "Deleting ModelContainer"
    ModelContainer.objects.all().delete()
    print "Deleting MetricContainer"
    MetricContainer.objects.all().delete()
    print "Deleting FeatureSelectionContainer"
    FeatureSelectionContainer.objects.all().delete()
    print "Deleting StatsModel"
    StatsModel.objects.all().delete()

    print "Deleting CatMolDescriptorValue"
    CatMolDescriptorValue.objects.exclude(
        descriptor__calculatorSoftware='manual').delete(recalculate_reactions=False)
    print "Deleting OrdMolDescriptorValue"
    OrdMolDescriptorValue.objects.exclude(
        descriptor__calculatorSoftware='manual').delete(recalculate_reactions=False)
    print "Deleting NumMolDescriptorValue"
    NumMolDescriptorValue.objects.exclude(
        descriptor__calculatorSoftware='manual').delete(recalculate_reactions=False)
    print "Deleting BoolMolDescriptorValue"
    BoolMolDescriptorValue.objects.exclude(
        descriptor__calculatorSoftware='manual').delete(recalculate_reactions=False)

    print "Deleting CatRxnDescriptorValue"
    CatRxnDescriptorValue.objects.exclude(
        descriptor__calculatorSoftware='manual').delete()

    print "Deleting NumRxnDescriptorValue"
    qs = NumRxnDescriptorValue.objects.exclude(
        descriptor__calculatorSoftware='manual')
    qs._raw_delete(qs.db)
    print "Deleting BoolRxnDescriptorValue"
    qs = BoolRxnDescriptorValue.objects.exclude(
        descriptor__calculatorSoftware='manual')
    qs._raw_delete(qs.db)
    print "Deleting OrdRxnDescriptorValue"
    qs = OrdRxnDescriptorValue.objects.exclude(
        descriptor__calculatorSoftware='manual')
    qs._raw_delete(qs.db)

    print "Deleting Descriptor"
    Descriptor.objects.exclude(calculatorSoftware='manual').delete()
