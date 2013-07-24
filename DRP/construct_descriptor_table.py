from models import *

#Take all data and make calculations.
def construct_entire_descriptor_table(lab_group):
	data_fields = get_model_field_names()
	for data in Data.objects.filter(lab_group=lab_group):
		compounds = [getattr(data,field) for field in data_fields if field[:-2]=="reactant"]
	
