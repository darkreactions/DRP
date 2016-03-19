import django
django.setup()
from DRP.models import PerformedReaction, DataSet, Descriptor, Compound, Reaction
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

#for d in Descriptor.objects.all():
    #if d.calculatorSoftware != "DRP":
        #print d.heading, d.calculatorSoftware

rxns = Reaction.objects.all()

rxns.calculate_descriptors(verbose=True)

#for i, r in enumerate(rxns):
    #if i > 1:
        #break
    #print "Calculating for reaction {}/{}: {}".format(i+1, rxns.count(), r)
    #r.save()

