#!/usr/bin/env python
import django
django.setup()
from DRP.models import ModelContainer, StatsModel

unbuilt_statsModel = StatsModel.objects.filter(container__built=False)

