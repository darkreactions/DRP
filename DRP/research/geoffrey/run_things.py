import django
django.setup()
from DRP.models import PerformedReaction, DataSet, Descriptor, Compound, Reaction, CompoundQuantity, NumRxnDescriptorValue
from django.db.models import Count

compounds = Compound.objects.all()

restart = 0

for i, c in enumerate(compounds):
    if i >= restart:
        print u"Calculating for compound {}/{}: {}".format(i+1, compounds.count(), c)
        try:
            c.save(calcDescriptors=True, invalidateReactions=False)
        except ValueError as e:
            print e

rxns = PerformedReaction.objects.all()

rxns = rxns.filter(valid=True).exclude(compounds=None)

rxns.calculate_descriptors(verbose=True)
