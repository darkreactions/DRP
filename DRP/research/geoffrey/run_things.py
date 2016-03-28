import django
django.setup()
from DRP.models import PerformedReaction, DataSet, BoolRxnDescriptorValue

rxns = PerformedReaction.objects.all()

print rxns.csvHeaders
