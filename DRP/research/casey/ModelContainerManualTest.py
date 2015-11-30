import os, sys
full_path = os.path.dirname(os.path.realpath(__file__))+"/"
django_path = full_path[:full_path.rfind("/DRP/")]
if django_path not in sys.path:
  sys.path = [django_path] + sys.path
  os.environ['DJANGO_SETTINGS_MODULE'] = 'DRP.settings'

from DRP.models import ModelContainer, Descriptor, PerformedReaction
from DRP.ml_models.splitters.KFoldSplitter import Splitter

#TODO: Get the predictors
responses = Descriptor.objects.filter(heading="crystallisation_outcome")
predictors = Descriptor.objects.exclude(id__in=responses.values_list('id', flat=True))

container = ModelContainer(library="weka", tool="svm", splitter="KFoldSplitter")
container.save()

reactions = PerformedReaction.objects.all()

splitter = Splitter()
for training, test in splitter.split(reactions):
  container.build(training, test, predictors, responses)

print container.summarize()


