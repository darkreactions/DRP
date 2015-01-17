from django.conf import settings
from data_config import CONFIG

from DRP.models import *

import json

#Import the data ranges from the json files.
static_dir = settings.BASE_DIR + settings.STATIC_URL
with open(static_dir+"js/editChoices.json") as file_handle:
 edit_choices = json.load(file_handle)
with open(static_dir+"js/dataRanges.json") as file_handle:
 data_ranges = json.load(file_handle)


#Create the form choices from the pre-defined ranges.
OUTCOME_CHOICES = [[opt,opt] for opt in edit_choices["outcomeChoices"]]
PURITY_CHOICES = [[opt,opt] for opt in edit_choices["purityChoices"]]
UNIT_CHOICES = [[opt,opt] for opt in edit_choices["unitChoices"]]
BOOL_CHOICES = [[opt,opt] for opt in edit_choices["boolChoices"]]
TYPE_CHOICES = [[opt,opt] for opt in edit_choices["typeChoices"]]

#Fields that are allowed to be stored as listy_strings.
list_fields = ["reactant", "quantity", "unit"]
#Fields that can be edited with a range alone: ###Copied in clientValidate.js
range_fields = {"quantity", "temp", "time", "pH", "outcome", "purity"}
#Fields that must be between a specific character count/limit:
limit_fields = {"ref", "notes"}
#Fields that must be a specific option:
opt_fields = {"unit", "slow_cool", "leak", "recommended"} #Note: slow_cool gains the general name of "slow"
bool_fields = {"slow_cool", "leak", "public", "is_valid", "recommended"}

#Type Groupings
int_fields = {"temp", "time",  "outcome", "purity"}
float_fields = {"quantity", "pH"}


def revalidate_datum(datum, lab_group):
 #Collect the data to validate
 dirty_data = {field:getattr(datum, field) for field in get_model_field_names()}
 #Validate and collect any errors
 (clean_data, errors) = full_validation(dirty_data, lab_group, revalidating=True)

 setattr(datum, "is_valid", clean_data["is_valid"])
 datum.save()

 return (clean_data, errors)

def full_validation(dirty_data, lab_group, revalidating=False):
 parsed_data = {} #Data that needs to be checked.
 clean_data = {} #Keep track of cleaned fields
 errors = {}

 fields = get_model_field_names()
 clean_data["is_valid"] = True

 #Gather the "coupled" fields (ie, the fields ending in a similar number)
 for field in list_fields:
  exec("{} = [[]]*{}".format(field, CONFIG.num_reactants)) #### {field: [[]]*CONFIG.num_reactants for field in list_fields} {field:
  parsed_data[field] = [[]]*CONFIG.num_reactants
  clean_data[field] = []
 ###fields = {field: [[]]*CONFIG.num_reactants for field in list_fields} ###CHANGE INTO ME, Future Casey
 ###parsed_data = {field: [[]]*CONFIG.num_reactants for field in list_fields}
 ###clean_data = {field: [] for field in list_fields}

 #Visible fields that are not required (not including rxn info).
 not_required = { ###Auto-generate?
  "notes", "duplicate_of"
 }

 for field in dirty_data:
  if field[-1].isdigit():
   #Put the data in its respective list.
   rel_list = eval("{}".format(field[:-2]))
   rel_list[int(field[-1])-1] = (dirty_data[field])
  else:
   try:
    assert(dirty_data[field]) #Assert that data was entered.
    parsed_data[field] = dirty_data[field]
   except:
    if field in not_required:
     clean_data[field] = "" #If nothing was entered, store nothing.
    else:
     errors[field] = "Field required."

 #Check that equal numbers of fields are present in each list
 for i in xrange(CONFIG.num_reactants):
  x = 0
  if reactant[i]:
   x+=2
   parsed_data["reactant"][i] = reactant[i]
  if quantity[i]:
   x+=3
   parsed_data["quantity"][i] = quantity[i]
  if x==5:
   parsed_data["unit"][i] = unit[i] #Menu, so no reason to check in form.

  #Unit is added automatically, so don't check it.
  if x == 3:
   errors["reactant_"+str(i+1)] = "Info missing."
  elif x == 2:
   errors["quantity_"+str(i+1)] = "Info missing."

 for field in parsed_data:
  #Make sure each reactant name is valid.
  if field=="reactant":
   for i in xrange(len(parsed_data[field])):
    if not parsed_data[field][i]: continue #Don't validate empty values.
    try:
     dirty_datum = str(parsed_data[field][i])
     if not clean_compound(dirty_datum)==dirty_datum:
      errors["{}_{}".format(field,i+1)] = "Contains illegal characters!"
     assert(validate_name(dirty_datum, lab_group))
     clean_data["{}_{}".format(field,i+1)] = dirty_datum #Add the filtered value to the clean values dict.
    except:
     errors["{}_{}".format(field,i+1)] = "Not in compound guide!"

  #Numeric fields:
  elif field in float_fields or field in int_fields:
   if field in float_fields: field_type="float"
   else: field_type="int"

   if field in list_fields:
    for i in xrange(len(parsed_data[field])):
     if not parsed_data[field][i]: continue #Don't validate empty values.
     try:
      dirty_datum = eval("{}(parsed_data[field][i])".format(field_type))
      assert(quick_validation(field, dirty_datum))
      clean_data["{}_{}".format(field,i+1)] = dirty_datum
     except:
      errors["{}_{}".format(field,i+1)] = "Must be between {} and {}.".format(data_ranges[field][0], data_ranges[field][1])
   else:
    try:
     dirty_datum = eval("{}(parsed_data[field])".format(field_type))
     assert(quick_validation(field, dirty_datum))
     parsed_data[field] = dirty_datum #Add the filtered mass to clean_data
     clean_data[field] = parsed_data[field]
    except:
     errors[field] = "Must be between {} and {}.".format(data_ranges[field][0], data_ranges[field][1])

  #Option fields:
  elif field in opt_fields:
   if field in list_fields:
    for i in xrange(len(parsed_data[field])):
     if not parsed_data[field][i]: continue #Don't validate empty values.
     try:
      dirty_datum = str(parsed_data[field][i])
      assert(quick_validation(field, dirty_datum))
      clean_data["{}_{}".format(field,i+1)] = dirty_datum
     except:
      if field in bool_fields:
       category="boolChoices"
      else:
       category = field+"Choices"

      errors["{}_{}".format(field,i+1)] = "Field must be one of: {}".format(edit_choices[category])
   else:
    try:
     dirty_datum = str(parsed_data[field])
     assert(quick_validation(field, dirty_datum))
     if clean_data["is_valid"]:
      clean_data["is_valid"] = CONFIG.unknown_label != dirty_datum
     clean_data[field] = dirty_datum
    except:
     if field in bool_fields:
      category="boolChoices"
     else:
      category = field+"Choices"

     errors[field] = "Field must be one of: {}".format(edit_choices[category])

  #Text fields.
  elif field in {"ref","notes", "duplicate_of"}:
   try:
    dirty_datum = str(parsed_data[field])
    assert(quick_validation(field, dirty_datum))

    #The "ref" already exists in the saved datum, so ignore it upon re-validation.
    if field=="ref" and not revalidating:
     try:
      #Gather the reference_set to make sure references are unique.
      ref_set = get_ref_set(lab_group)

      #Check to make sure the ref isn't in the ref_set.
      assert(not dirty_datum in ref_set)
      clean_data[field] = dirty_datum
     except:
      errors[field] = "Already in use."

    elif field=="duplicate_of":
     try:
      #Gather the reference_set to make sure references are unique.
      ref_set = get_ref_set(lab_group)

      #Check to make sure the ref isn't in the ref_set.
      assert(dirty_datum in ref_set)
      clean_data[field] = dirty_datum
     except:
      errors[field] = "Nonexistent reference."
    else:
     clean_data[field] = dirty_datum
   except:
    errors[field] = "Cannot exceed {} characters.".format(data_ranges[field][1])

 return (clean_data, errors)


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


def validate_name(abbrev_to_check, lab_group):
 #Get the cached set of abbreviations.
 abbrevs = collect_CG_name_pairs(lab_group)
 return abbrev_to_check in abbrevs


def clean_compound(compound):
  #Remove any non-printable characters.
  return filter(lambda x: x in string.printable, compound)


def validate_CG(dirty_data, lab_group, editing_this=False):
  def parse_CAS_ID(CAS):
    CAS = CAS.replace(" ", "-").replace("/", "-").replace("_", "-")

    #Check that the CAS ID has three hyphen-delineated parts.
    if len(CAS.split("-")) != 3:
      raise Exception("CAS ID requires three distinct parts.")

    #Check that only numbers are present.
    elif not CAS.replace("-","").isdigit():
      raise Exception("CAS ID may only have numeric characters.")

    return CAS

  from DRP.chemspider import search_chemspider

  #Variable Setup
  clean_data = dirty_data
  errors = {}

  for field in ["compound", "abbrev", "compound_type"]:
    if not dirty_data.get(field):
      errors[field] = "Field required."

  #Get the CAS_ID if applicable.
  raw_CAS = dirty_data.get("CAS_ID")
  try:
    clean_data["CAS_ID"] = parse_CAS_ID(raw_CAS) if raw_CAS else ""
  except Exception as e:
    errors["CAS_ID"] = e

  #If the data is custom, don't query ChemSpider.
  if dirty_data.get("custom"):
    clean_data["custom"]=True
    clean_data["image_url"], clean_data["smiles"], clean_data["mw"] = "","",""
  else: #But if it is normal, get extra data from the query.
    try:
      search_fields = [dirty_data.get("CAS_ID"), dirty_data.get("compound")]
      query = search_chemspider(search_fields)
      clean_data["image_url"] = query.imageurl
      clean_data["smiles"] = query.smiles
      clean_data["mw"] = query.molecularweight

    except:
      if search_fields[0]:
        errors["CAS_ID"] = "Could not find a molecule with this CAS ID."
      else:
        errors["compound"] = "Could not find this compound."

  #Prevent duplicate abbrevs.

  if not errors.get("abbrev"):
    clean_data["abbrev"] = dirty_data["abbrev"]
    if not editing_this:
      if get_lab_CG(lab_group).filter(abbrev=clean_data["abbrev"]).exists():
        errors["abbrev"] = "Abbreviation already used."

  return clean_data, errors
