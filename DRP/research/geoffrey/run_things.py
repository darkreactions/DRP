#!/usr/bin/env python
import django
django.setup()
import DRP
from itertools import chain

headings = []

for compoundRole in DRP.models.CompoundRole.objects.all():
    for descriptor in chain(DRP.models.BoolMolDescriptor.objects.filter(heading__contains='group'), DRP.models.BoolMolDescriptor.objects.filter(heading__contains='period'), DRP.models.BoolMolDescriptor.objects.filter(heading__contains='valence')):
        for i in (True, False):
            heading = '{}_{}_{}_count'.format(compoundRole.label, descriptor.csvHeader, i)
            headings.append(heading)
            heading = '{}_{}_{}_molarity'.format(compoundRole.label, descriptor.csvHeader, i)
            headings.append(heading)
        heading = '{}_{}_any'.format(compoundRole.label, descriptor.csvHeader)
        headings.append(heading)

print '\n'.join(headings)
    
