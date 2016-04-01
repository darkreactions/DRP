import django
django.setup()
from DRP.models import ModelContainer
import build_model
from sys import argv

model_name = argv[1]
desc_tail = ' BCR Weighted. Legacy rxns. New descriptors with CA at rxn pH NonzeroVariance. 15 mutual info split.'
desc = model_name + desc_tail

m = ModelContainer.objects.get(description=desc)

print m.description
print 'pk: {}'.format(m.pk)
print m.splitter
print m.modelVisitorLibrary, m.modelVisitorTool

build_model.display_model_results(m)
