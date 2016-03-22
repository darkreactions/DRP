import django
django.setup()
from DRP.models import PerformedReaction, DataSet, Descriptor, Compound, Reaction, CompoundQuantity, NumRxnDescriptorValue, CatMolDescriptorValue, OrdMolDescriptorValue, BoolMolDescriptorValue, NumMolDescriptorValue, CatMolDescriptor, OrdMolDescriptor, BoolMolDescriptor, NumMolDescriptor, BoolRxnDescriptorValue
from django.db.models import Count

import argparse

parser = argparse.ArgumentParser(description='dumps to csv', fromfile_prefix_chars='@',
                                 epilog="Prefix arguments with '@' to specify a file containing newline"
                                 "-separated values for that argument. e.g.'-p @predictor_headers.txt'"
                                 " to pass multiple descriptors from a file as predictors")

parser.add_argument('-p', '--predictor-headers', nargs='+',
                    help='One or more descriptors to use as predictors.', required=True)
args = parser.parse_args()

rxns = DataSet.objects.get(name='valid_legacy_reactions').reactions.all()
print rxns.count()
with open('reactions_old.csv', 'wb') as f:
    rxns.toCsv(f, expanded=True, whitelistHeaders=args.predictor_headers)
