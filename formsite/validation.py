###GAE serves static files separately and thus they can't easily be accessed.
##from django.conf import settings
##import json
##with open(settings.STATIC_ROOT+"js/editChoices.json") as file_handle:
	##range_fields = json.load(file_handle)

#Contents of JSON files:
data_ranges = {
	"quantity" : [0,50],
	"temp" : [0,500],
	"outcome" : [0,4],
	"purity" : [0,2],
	"time" : [0,350],
	"pH" : [-2,17],
	"ref" : [1,12],
	"notes" : [0,65]
	}
edit_choices = {
	"unitChoices" : ["g","mL","d"],
	"boolChoices" : ["Yes","No","?"],
	"outcomeChoices" : [0,1,2,3,4],
	"purityChoices" : [0,1,2],
	"typeChoices": ["Org", "Inorg", "pH", "Ox", "Sol", "Water"],
	}
	
#Fields that can be edited with a range alone: ###Copied in clientValidate.js
range_fields = {"quantity", "temp", "time", "pH", "outcome", "purity"}
#Fields that must be between a specific character count/limit:
limit_fields = {"ref", "notes"}
#Fields that must be a specific option:
opt_fields = {"unit", "slow_cool", "leak"} #Note: slow_cool gains the general name of "slow" 
bool_fields = {"slow_cool", "leak"}

#Type Groupings
int_fields = {"temp", "time",  "outcome", "purity"}
float_fields = {"quantity", "pH"}

			
def quick_validation(field, dirty_datum, model="Data"):
	if model=="Data":
		if field in data_ranges:
			data_range = data_ranges[field]
			
		if field in range_fields:
			return (data_range[0] <= float(dirty_datum) <= data_range[1])
		if field in limit_fields:
			return (data_range[0] <= len(dirty_datum) <= data_range[1])
		if field in opt_fields:
			if field in bool_fields: category = "boolChoices"
			else: category = field+"Choices"
			return dirty_datum in edit_choices[category]
		return True
	if model=="CompoundEntry":
		if field=="CAS_ID":
			length_valid = (len(dirty_datum.split("-"))==3)
			elements_valid = dirty_datum.translate(None, "-").isdigit()
			return length_valid and elements_valid
		if field=="compound_type":
			return dirty_datum in edit_choices["typeChoices"] 
		return True
