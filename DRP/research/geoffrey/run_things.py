#!/usr/bin/env python
import django
django.setup()
import DRP
from DRP.models import Descriptor, NumRxnDescriptor, NumMolDescriptor
from itertools import chain


headings = [d.heading for d in NumRxnDescriptor.objects.filter(
    heading__contains='pHreaction')]

print '\n'.join(headings)
