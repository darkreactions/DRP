#!/usr/bin/env python
import django
django.setup()
from DRP.models import ModelContainer
import build_model
from sys import argv


m = ModelContainer.objects.order_by('-pk')[0]

new_m = m.create_duplicate()

new_m.build(verbose=True)
