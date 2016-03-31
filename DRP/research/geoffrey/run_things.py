import django
django.setup()
from DRP.models import Descriptor, Compound


with open('test.arff', 'w') as f:
    Compound.objects.all().toArff(f)
