#!/usr/bin/env python
import django
django.setup()
import DRP
from DRP.models import Descriptor


print '\n'.join([d.heading for d in Descriptor.objects.filter(heading__contains='polar_surface_area')])
