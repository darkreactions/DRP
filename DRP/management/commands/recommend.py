from django.core.management.base import BaseCommand
from DRP.recommender.ReactionRecommender import Recommender
from DRP.models import LabGroup, Reaction, PerformedReaction, Descriptor, ModelContainer, NumRxnDescriptorValue, NumRxnDescriptor, BoolRxnDescriptor
from DRP.plugins.rxndescriptors.rxnhash import calculate_many
from DRP.ml_models.splitters.kFoldSplitter import Splitter
from django.core.exceptions import ObjectDoesNotExist

import numpy as np

def get_numeric_descriptor_amounts(num_desc_ids, labGroup_id):
	""" Take a list of reaction descriptor ids and return a dictionary of descriptor ids to reasonable amounts based
		off of reactions performed by a the lab specified by the labGroup_id"""
	reasonable_desc_values_dict = {}
	performed_reaction_ids = PerformedReaction.objects.all().values_list('id', flat=True).order_by('id')
	reaction_ids = Reaction.objects.filter(labGroup_id=labGroup_id, id__in=performed_reaction_ids).values_list('id', flat=True).order_by('id')
	print(len(reaction_ids))
	for num_desc_id in num_desc_ids:
		values = NumRxnDescriptorValue.objects.filter(descriptor_id=num_desc_id, reaction_id__in=reaction_ids).exclude(value=None).values_list('value', flat=True)
		print("Values ", values)
		values = np.array(values)
		mean_value = np.mean(values)
		value_sd = np.std(values)
		reasonable_values = [mean_value - 2 * value_sd, mean_value - value_sd, mean_value, mean_value + value_sd, mean_value + 2 * value_sd]
		reasonable_values = [value for value in reasonable_values if value > 0]
		reasonable_desc_values_dict[num_desc_id] = reasonable_values
	return reasonable_desc_values_dict



class Command(BaseCommand):

	help = 'Recommends a set of reactions to be performed'

	def add_arguments(self, parser):
		parser.add_argument('labGroup', type=str)
	
	def handle(self, *args, **options):

		# TODO: make labGroup_id an argument to pass in based off of labGroup's name
		try:
			labGroup_id = LabGroup.objects.get(title=options['labGroup'])
		except ObjectDoesNotExist:
			print('Error: "{}" is not a valid lab group.'.format(options['labGroup']))
			exit()


		allNumericDescriptors = NumRxnDescriptor.objects.all()
		boolean_crystallisation_outcome_id = Descriptor.objects.get(heading='boolean_crystallisation_outcome', calculatorSoftware='manual').id
		boolean_crystallisation_outcome_descriptor_object = BoolRxnDescriptor.objects.get(booleandescriptor_ptr_id=boolean_crystallisation_outcome_id)
		
		desired_desc_dict = { boolean_crystallisation_outcome_descriptor_object : [1]}

		''' descriptors for grid_params needed: 
			- reaction_temperature
			- reaction_time
			- reaction_pH
			- reaction pre heat standing time
			- was this reaction performed in a teflon pouch?
			- did this reaction leak?
			- was a slow cool performed?
			- was this reaction performed in an oil bath?
		'''

		# necessary numeric reactions
		reaction_temperature_id = 	Descriptor.objects.get(heading='reaction_temperature').id
		reaction_time_id = 			Descriptor.objects.get(heading='reaction_time').id
		reaction_pH_id = 			Descriptor.objects.get(heading='reaction_pH').id
		pre_heat_standing_id = 		Descriptor.objects.get(heading='pre_heat_standing').id

		reaction_temperature =		NumRxnDescriptor.objects.get(numericdescriptor_ptr_id=reaction_temperature_id)
		reaction_time =				NumRxnDescriptor.objects.get(numericdescriptor_ptr_id=reaction_time_id)
		reaction_pH =				NumRxnDescriptor.objects.get(numericdescriptor_ptr_id=reaction_pH_id)
		pre_heat_standing =			NumRxnDescriptor.objects.get(numericdescriptor_ptr_id=pre_heat_standing_id)

		teflon_pouch_id = 			Descriptor.objects.get(heading='teflon_pouch').id
		leak_id = 					Descriptor.objects.get(heading='leak').id
		slow_cool_id = 				Descriptor.objects.get(heading='slow_cool').id
		oil_bath_id = 				Descriptor.objects.get(heading='oil_bath').id

		teflon_pouch =				BoolRxnDescriptor.objects.get(booleandescriptor_ptr_id=teflon_pouch_id)
		leak =						BoolRxnDescriptor.objects.get(booleandescriptor_ptr_id=leak_id)
		slow_cool =					BoolRxnDescriptor.objects.get(booleandescriptor_ptr_id=slow_cool_id)
		oil_bath =					BoolRxnDescriptor.objects.get(booleandescriptor_ptr_id=oil_bath_id)

		numeric_descriptor_amounts = get_numeric_descriptor_amounts([reaction_temperature_id, reaction_time_id, reaction_pH_id, pre_heat_standing_id], labGroup_id)

		grid_params = {
						reaction_temperature: numeric_descriptor_amounts[reaction_temperature_id],
						reaction_time: numeric_descriptor_amounts[reaction_time_id],
						reaction_pH: numeric_descriptor_amounts[reaction_pH_id],
						pre_heat_standing: numeric_descriptor_amounts[pre_heat_standing_id],
						teflon_pouch: [0, 1],
						leak: [0, 1],
						slow_cool: [0, 1],
						oil_bath: [0, 1] }



		print("Initialize recommender")

		# TODO: address magic number from this next line
		obj = ModelContainer.objects.filter(modelVisitorTool = "SVM_PUK")[4]
		print(obj.built)

		rec = Recommender(obj, grid_params, desired_desc_dict)

		print("Starting recommendation")
		reactions = rec.recommend()
		print(reactions)
