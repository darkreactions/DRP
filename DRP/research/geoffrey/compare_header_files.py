import django
django.setup()
from DRP.models import PerformedReaction, DataSet, Descriptor, Compound, Reaction, CompoundQuantity, NumRxnDescriptorValue, CatMolDescriptorValue, OrdMolDescriptorValue, BoolMolDescriptorValue, NumMolDescriptorValue, CatMolDescriptor, OrdMolDescriptor, BoolMolDescriptor, NumMolDescriptor, BoolRxnDescriptorValue
from django.db.models import Count

import argparse

parser = argparse.ArgumentParser(description='dumps to csv', fromfile_prefix_chars='@',
                                 epilog="Prefix arguments with '@' to specify a file containing newline"
                                 "-separated values for that argument. e.g.'-p @predictor_headers.txt'"
                                 " to pass multiple descriptors from a file as predictors")

parser.add_argument('-o', '--headers1', nargs='+', required=True)
parser.add_argument('-n', '--headers2', nargs='+', required=True)
args = parser.parse_args()

missing_headers = [h for h in args.headers1 if h not in args.headers2]
for header in missing_headers:
    print header
