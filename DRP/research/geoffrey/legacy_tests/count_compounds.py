import django
django.setup()
from DRP.models import PerformedReaction, DataSet
from django.db.models import Count

reactions = PerformedReaction.objects.annotate(compound_count=Count('compounds'))
reactions = reactions.filter(valid=True)
compound_numbers = [reactions.filter(compound_count=i).count() for i in xrange(10)]
print compound_numbers
