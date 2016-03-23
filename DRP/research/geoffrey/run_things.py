import django
django.setup()
from DRP.models import PerformedReaction, DataSet, Descriptor, Compound, Reaction, CompoundQuantity, NumRxnDescriptorValue, CatMolDescriptorValue, OrdMolDescriptorValue, BoolMolDescriptorValue, NumMolDescriptorValue, CatMolDescriptor, OrdMolDescriptor, BoolMolDescriptor, NumMolDescriptor, BoolRxnDescriptorValue
from django.db.models import Count

# compounds = Compound.objects.all()

# compounds.calculate_descriptors(verbose=True)

rxns = PerformedReaction.objects.all()

rxns.calculate_descriptors(verbose=True)
