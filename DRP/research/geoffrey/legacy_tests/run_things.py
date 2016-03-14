import django
django.setup()
from DRP.models import PerformedReaction, DataSet, Descriptor
from django.db.models import Count

for d in Descriptor.objects.all():
    print d.heading
