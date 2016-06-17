import os
import sys
full_path = os.path.dirname(os.path.realpath(__file__)) + "/"
django_path = full_path[:full_path.rfind("/DRP/")]
if django_path not in sys.path:
    sys.path = [django_path] + sys.path
    os.environ['DJANGO_SETTINGS_MODULE'] = 'DRP.settings'

from DRP.models import ModelContainer, Descriptor, PerformedReaction
from DRP.ml_models.splitters.KFoldSplitter import Splitter

# TODO: Get the predictors
responses = Descriptor.objects.filter(heading="crystallisation_outcome")
predictors = Descriptor.objects.exclude(
    id__in=responses.values_list('id', flat=True))[:15]

container = ModelContainer(
    library="weka", tool="svm", splitter="KFoldSplitter")
container.save()

#reactions = PerformedReaction.objects.all()
reactions = PerformedReaction.objects.filter(reference__istartswith="1")

print "# of Reactions: {}".format(reactions.count())
print "# of Predictors: {}".format(predictors.count())
print "# of Responses: {}".format(responses.count())

splitter = Splitter()
for training, test in splitter.split(reactions):
    print "\nBuilding a model!"
    container.build(training, test, predictors, responses, debug=True)

print container.summarize()
