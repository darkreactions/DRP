import django
django.setup()
from DRP.models import Descriptor
from django.core.exceptions import ValidationError


calcSoftMap = [
    ('example.py plugin', 'example_plugin'),
    ('ChemAxon/cxcalc', 'ChemAxon_cxcalc'),
    ('drp/rdkit', 'drp_rdkit'),
    ('DRP/xxhash', 'DRP_xxhash'),
]

calcSoftVersMap = [
    ('0.02/0.5.0', '0.02_0.5.0'),
    ('0.02/0.4.3', '0.02_0.4.3')
]

print Descriptor.objects.values('calculatorSoftware').distinct()

for old, new in calcSoftMap:
    Descriptor.objects.filter(calculatorSoftware=old).update(
        calculatorSoftware=new)

print Descriptor.objects.values('calculatorSoftware').distinct()

print Descriptor.objects.values('calculatorSoftwareVersion').distinct()

for old, new in calcSoftVersMap:
    Descriptor.objects.filter(calculatorSoftwareVersion=old).update(
        calculatorSoftwareVersion=new)

print Descriptor.objects.values('calculatorSoftwareVersion').distinct()

for d in Descriptor.objects.all():
    try:
        d.full_clean()
    except ValidationError:
        print d.name, d.heading, d.calculatorSoftware, d.calculatorSoftwareVersion
        raise
