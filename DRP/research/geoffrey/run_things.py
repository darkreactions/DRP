#!/usr/bin/env python
import django
django.setup()
import DRP
from DRP.models import Descriptor, NumRxnDescriptor, NumMolDescriptor
from itertools import chain


headings = chain([d.heading for d in NumMolDescriptor.objects.filter(heading__contains='nominal')], [d.heading for d in NumRxnDescriptor.objects.filter(heading__contains='nominal')])

print '\n'.join(headings)
