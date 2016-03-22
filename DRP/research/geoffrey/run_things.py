import django
django.setup()
from DRP.models import PerformedReaction, DataSet, Descriptor, Compound, Reaction, CompoundQuantity, NumRxnDescriptorValue, CatMolDescriptorValue, OrdMolDescriptorValue, BoolMolDescriptorValue, NumMolDescriptorValue, CatMolDescriptor, OrdMolDescriptor, BoolMolDescriptor, NumMolDescriptor, BoolRxnDescriptorValue
from django.db.models import Count

rxns = PerformedReaction.objects.all()
with open('reactions_old.csv', 'wb') as f:
    rxns.toCsv(f, expanded=True)
