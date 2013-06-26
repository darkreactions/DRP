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
	"ref" : [0,12],
	"notes" : [0,65]
	}
edit_choices = {
	"unitChoices" : ["g","mL","d"],
	"boolChoices" : ["Yes","No","?"],
	"outcomeChoices" : [0,1,2,3,4],
	"purityChoices" : [0,1,2]
	}
	
#Fields that can be edited with a range alone: ###Copied in clientValidate.js
range_fields = {"quantity", "temp", "time", "pH", "outcome", "purity"}
#Fields that must be between a specific character count/limit:
limit_fields = {"ref", "notes"}
def quick_validation(field, dirty_datum):
	if field in data_ranges:
		data_range = data_ranges[field]
		
	if field in range_fields:
		return (data_range[0] <= float(dirty_datum) <= data_range[1])
	elif field in limit_fields:
		return (data_range[0] <= len(dirty_datum) <= data_range[1])
	else:
		return True
	
