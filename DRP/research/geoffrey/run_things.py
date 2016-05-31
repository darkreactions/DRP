#!/usr/bin/env python
import django
django.setup()
from DRP.models import *
from itertools import chain


for d in BoolRxnDescriptor.objects.filter(heading__contains='valence'):
    print d.heading
    d.heading = d.heading.replace('boolean_inorganic_valence', 'drpInorgAtom_boolean_valence')
    d.save()
