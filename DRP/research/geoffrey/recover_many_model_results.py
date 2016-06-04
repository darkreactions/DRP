import django
django.setup()
from DRP.models import ModelContainer
import build_model
from sys import argv

model_names = argv[1:]
desc = "Legacy rxns. CFS of new descriptors with ChemAxon at reaction pH. 15 mutual info split."

containers = ModelContainer.objects.filter(description__startswith=desc, built=True)
print containers.count()

for m in containers:
    print m.description
    print 'pk: {}'.format(m.pk)
    print m.splitter
    print m.modelVisitorLibrary, m.modelVisitorTool

    build_model.display_model_results(m, heading=m.modelVisitorTool + m.modelVisitorOptions)
