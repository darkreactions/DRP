import django
django.setup()
from DRP.models import (DataSet, DataSetRelation, ModelContainer, MetricContainer, FeatureSelectionContainer, StatsModel,
                     CatRxnDescriptorValue, OrdRxnDescriptorValue, BoolRxnDescriptorValue, NumRxnDescriptorValue, 
                     CatMolDescriptorValue, OrdMolDescriptorValue, BoolMolDescriptorValue, NumMolDescriptorValue, 
                     CatRxnDescriptor, OrdRxnDescriptor, BoolRxnDescriptor, NumRxnDescriptor, 
                     CatMolDescriptor, OrdMolDescriptor, BoolMolDescriptor, NumMolDescriptor,
                     RxnDescriptorValue,
                     Descriptor)
                    

models_to_delete = [DataSet, DataSetRelation, ModelContainer, MetricContainer, FeatureSelectionContainer, StatsModel,
                     CatRxnDescriptorValue, OrdRxnDescriptorValue, BoolRxnDescriptorValue, NumRxnDescriptorValue, 
                     CatMolDescriptorValue, OrdMolDescriptorValue, BoolMolDescriptorValue, NumMolDescriptorValue, 
                     CatRxnDescriptor, OrdRxnDescriptor, BoolRxnDescriptor, NumRxnDescriptor, 
                     CatMolDescriptor, OrdMolDescriptor, BoolMolDescriptor, NumMolDescriptor,
                     Descriptor
                    ]

special_instruction_models = [CatMolDescriptorValue, OrdMolDescriptorValue, BoolMolDescriptorValue, NumMolDescriptorValue]
                    
if __name__ == '__main__':
    really = raw_input("Are you sure you want to delete all but the raw data in this database? This includes all models, descriptors, descriptor values and containers and is irreversible: ")
    if really.lower() == 'yes':
        for m in models_to_delete:
            qs = m.objects.all()
            print "Deleting {} {} objects".format(qs.count(), m)
            if m in special_instruction_models:
                qs.delete(recalculate_reactions=False)
            elif qs.count() > 500000: # We can only delete so many things at once. Filtering is better than slicing for speed too, I think
                pk_cutoff = 500000
                while qs.count() > 0:
                    print "Deleting up to pk {}".format(pk_cutoff)
                    qs.filter(pk__lte=pk_cutoff).delete()
                    qs = m.objects.all()
                    pk_cutoff += 500000
            else:
                qs.delete()
        # sweep again
        for m in models_to_delete:
            qs = m.objects.all()
            if qs.count() != 0:
                raise RuntimeError("Failed to delete all objects of model {}. {} remain".format(m, qs.count()))
