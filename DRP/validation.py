from django.conf import settings
from data_config import CONFIG
import json

#Import the data ranges from the json files.
static_dir = settings.BASE_DIR + settings.STATIC_URL
with open(static_dir+"js/editChoices.json") as file_handle:
 edit_choices = json.load(file_handle)
with open(static_dir+"js/dataRanges.json") as file_handle:
 data_ranges = json.load(file_handle)


#Fields that can be edited with a range alone: ###Copied in clientValidate.js
range_fields = {"quantity", "temp", "time", "pH", "outcome", "purity"}
#Fields that must be between a specific character count/limit:
limit_fields = {"ref", "notes"}
#Fields that must be a specific option:
opt_fields = {"unit", "slow_cool", "leak", "recommended"} #Note: slow_cool gains the general name of "slow"
bool_fields = {"slow_cool", "leak", "public", "valid", "recommended"}

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
