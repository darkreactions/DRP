from django.core.management.base import BaseCommand
from DRP.recommender.ReactionRecommender import Recommender
from DRP.models import PerformedReaction, Descriptor, ModelContainer, NumRxnDescriptorValue, NumRxnDescriptor
from DRP.plugins.rxndescriptors.rxnhash import calculate_many
from DRP.ml_models.splitters.kFoldSplitter import Splitter

class Command(BaseCommand):

	help = 'Reccomends a set of reactions to be performed'

	def handle(self, *args, **options):
		
		# responses = Descriptor.objects.filter(heading="crystallisation_outcome")
		# predictors = Descriptor.objects.exclude(id__in=responses.values_list('id', flat=True))[:20]
		# reactions = PerformedReaction.objects.all()
		# splitter = Splitter("kFoldSplitter")
		# container = ModelContainer.create(modelVisitorLibrary="weka", modelVisitorTool="SVM_PUK", \
		# 	reactions=reactions, predictors=predictors, responses=responses, splitter="kFoldSplitter")
		
		# print("Created container successfully")

		# container.build(verbose = True)

		# print("Built container successfully")

		allNumericDescriptors = NumRxnDescriptor.objects.all()

		# values = NumRxnDescriptorValue.objects.all()
		# values = values[:2]
		# values = [val.value for val in values]


		grid_params = {allNumericDescriptors[0] : [0, 1, 2, 3, 4, 5]}

		desired_desc_dict = { allNumericDescriptors[0] : [0, 1, 2, 3, 4, 5]}

		print("Initialize reccomender")

		obj = ModelContainer.objects.filter(modelVisitorTool = "SVM_PUK")[0]
		print(obj.built)
		# obj.build()
		# obj.save()

		rec = Recommender(obj, grid_params, desired_desc_dict)

		print("Starting recommendation")
		reactions = rec.recommend()
		print(reactions)
