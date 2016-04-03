#!/usr/bin/env python
import django
django.setup()
from DRP.models import ModelContainer
import build_model
from sys import argv


m = ModelContainer.objects.filter(statsmodel__isnull=False).order_by('-pk')[0]

options = {'BCR': False}

new_m = m.create_duplicate(modelVisitorTool='J48', modelVisitorOptions=options)

new_m.build(verbose=True)
