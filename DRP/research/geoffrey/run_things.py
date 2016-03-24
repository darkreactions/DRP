import django
django.setup()
from DRP.models import PerformedReaction, DataSet, Descriptor, Compound, Reaction, CompoundQuantity, NumRxnDescriptorValue, CatMolDescriptorValue, OrdMolDescriptorValue, BoolMolDescriptorValue, NumMolDescriptorValue, CatMolDescriptor, OrdMolDescriptor, BoolMolDescriptor, NumMolDescriptor, BoolRxnDescriptorValue
from django.db.models import Count
import argparse

#print PerformedReaction.objects.all().expandedArffHeaders

print DataSet.object.get(name='test1001').reactions.all().expandedArffHeaders

#if __name__ == '__main__':
    #django.setup()
    #parser = argparse.ArgumentParser(description='Builds a model', fromfile_prefix_chars='@',
                                     #epilog="Prefix arguments with '@' to specify a file containing newline"
                                     #"-separated values for that argument. e.g.'-p @predictor_headers.txt'"
                                     #" to pass multiple descriptors from a file as predictors")
    #parser.add_argument('-d', '--descriptor-headers', nargs='+',
                        #help='One or more descriptors to use as predictors.', required=True)
    #parser.add_argument('-rxn', '--reaction-set-name', default=None,
                        #help='The name of the reactions to use as a whole dataset')
    #args = parser.parse_args()

    #print args.descriptor_headers
