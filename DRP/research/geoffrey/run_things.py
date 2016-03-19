import django
django.setup()
from DRP.models import PerformedReaction, DataSet, Descriptor, Compound, Reaction, CompoundQuantity, NumRxnDescriptorValue
from django.db.models import Count

compounds = Compound.objects.all()
print compounds.count()

restart = 0

compounds.calculate_descriptors(verbose=True)

#bad_compounds = []

#for i, c in enumerate(compounds):
    #if i+1 >= restart:
        #print u"Calculating for compound {}/{}: {}".format(i+1, compounds.count(), c)
        ##try:
        #c.save(calcDescriptors=True, invalidateReactions=False)
        ##except ValueError as e:
            ##bad_compounds.append(c)
            ##print e

#print "These compounds were not found: "
#for c in bad_compounds:
    #try:
        #print c
    #except UnicodeEncodeError as e:
        #print "Error on compound with abbreviation: ", c.abbrev
        #print e

#rxns = PerformedReaction.objects.all()
#rxns = rxns.filter(valid=True).exclude(compounds=None)
#rxns.calculate_descriptors(verbose=True)

