#!/usr/bin/env python
import django
django.setup()
from DRP.models import Descriptor

for d in Descriptor.objects.filter(heading__contains='valence'):
    print d.heading
