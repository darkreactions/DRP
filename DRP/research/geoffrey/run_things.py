import django
django.setup()
from DRP.models import PerformedReaction, DataSet, Descriptor, Compound, Reaction, CompoundQuantity, NumRxnDescriptorValue, CatMolDescriptorValue, OrdMolDescriptorValue, BoolMolDescriptorValue, NumMolDescriptorValue, CatMolDescriptor, OrdMolDescriptor, BoolMolDescriptor, NumMolDescriptor
from django.db.models import Count

compounds = Compound.objects.all()
###restart = 0

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


#print NumRxnDescriptorValue.objects.exclude(value=None).count()

#for cmd in CatMolDescriptor.objects.all().order_by('heading'):
    #all_vals = CatMolDescriptorValue.objects.filter(descriptor=cmd)
    #print cmd, cmd.calculatorSoftware, all_vals.count(), all_vals.exclude(value=None).count()
#for nmd in NumMolDescriptor.objects.all().order_by('heading'):
    #all_vals = NumMolDescriptorValue.objects.filter(descriptor=nmd)
    #print nmd, nmd.calculatorSoftware, all_vals.count(), all_vals.exclude(value=None).count()
#for omd in OrdMolDescriptor.objects.all().order_by('heading'):
    #print omd, omd.calculatorSoftware, OrdMolDescriptorValue.objects.filter(descriptor=omd).exclude(value=None).count()
#for bmd in BoolMolDescriptor.objects.all().order_by('heading'):
    #print bmd, bmd.calculatorSoftware, BoolMolDescriptorValue.objects.filter(descriptor=bmd).exclude(value=None).count()
