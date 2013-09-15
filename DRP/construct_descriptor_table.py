from models import *
from data_config import CONFIG

#Take ALL data and make calculations.
def construct_entire_descriptor_table(lab_group):
	#Variable Setup.
	skip = {}
	CG_translate = {}
	abbrev_to_type = {}
	to_skip = {}
	
	#Gather the CG entries.
	for CG in CompoundEntry.objects.filter(lab_group=lab_group):
		try:
			CG_translate[CG.abbrev] = CG.compound
			abbrev_to_type[CG.abbrev] = CG.compound_type
			
		except:
			print "Error: {} failed to load!".format(CG.compound)
	
	#Gather the valid Data entries that don't have calculations.
	data_fields = get_model_field_names()
	count = 0
	for data in Data.objects.filter(lab_group=lab_group, calculations=None, is_valid=True):
		#Display counter information:
		if count % 100000 == 0:
			print "---- Played with {} reactions so far...".format(count)
		count += 1		
		
		#Collect the reactant information.
		compounds = [getattr(data,field) for field in data_fields if field[:-2]=="reactant"]
		quantities = [getattr(data,field) for field in data_fields if field[:-2]=="quantity"]
		units = [getattr(data,field) for field in data_fields if field[:-2]=="unit"]
		
		for unit in units:
			if unit != "g":
				###Convert any amounts that need to be converted.
				pass###
		
		###Assume data is valid?
		if not data.is_valid:
			pass###
			
		
		
	
