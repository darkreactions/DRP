from models import *
from data_config import CONFIG

#Take ALL data and make calculations.
def construct_entire_descriptor_table(lab_group):
	#Variable Setup.
	skip = {}
	CG_translate = {}
	
	#Gather the CG entries.
	for CG in CompoundEntry.objects.filter(lab_group=lab_group):
		CG_translate[CG.abbrev] = CG.compound
	
	#Gather the Data entries that don't have calculations.
	data_fields = get_model_field_names()
	for data in Data.objects.filter(lab_group=lab_group, calculations=None):
		#Collect the reactant information.
		compounds = [getattr(data,field) for field in data_fields if field[:-2]=="reactant"]
		quantities = [getattr(data,field) for field in data_fields if field[:-2]=="quantity"]
		units = [getattr(data,field) for field in data_fields if field[:-2]=="unit"]
		
		for unit in units:
			if unit != "g":
				###Convert any amounts that need to be converted.
				pass###
		
		
	
