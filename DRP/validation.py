from django.conf import settings
from data_config import CONFIG

from DRP.models import get_model_field_names

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
allow_unknown = {"outcome","purity","slow_cool","leak","recommended"}

#Type Groupings
int_fields = {"temp", "time"}
float_fields = {"quantity", "pH"}


def revalidate_datum(datum, lab_group):
 #Collect the data to validate
 dirty_data = {field:getattr(datum, field) for field in get_model_field_names()}
 #Validate and collect any errors
 (clean_data, errors) = full_validation(dirty_data, lab_group, revalidating=True)

 setattr(datum, "is_valid", clean_data["is_valid"])
 datum.save()

 return (clean_data, errors)


def is_numeric(val):
  try:
    float(val)
    return True
  except:
    return False



def full_validation(dirty_data, lab_group, revalidating=False):
  from DRP.models import compound_exists, get_compound

  clean_data = {} #Keep track of cleaned fields
  errors = {}

  clean_data["is_valid"] = True

  for i in CONFIG.reactant_range():
    clean_data["unit_{}".format(i)] = "g"

  reactants = []
  used_abbrevs = set()

  for i in CONFIG.reactant_range():
    reactant = "reactant_fk_{}".format(i)
    if reactant in dirty_data:
      abbrev = dirty_data[reactant]

      if abbrev in used_abbrevs:
        errors[reactant] = "Compound already used!"

      if not compound_exists(abbrev, lab_group=lab_group):
        errors[reactant] = "Compound not found!"

      if "quantity_{}".format(i) in dirty_data:
        quantity = dirty_data["quantity_{}".format(i)]

        if not is_numeric(quantity):
          errors["quantity_{}".format(i)] = "Must be a numeric value."

        if "unit_{}".format(i) in dirty_data:
          unit = dirty_data["unit_{}".format(i)]
          used_abbrevs.add(abbrev)

          reactants.append({
            "reactant":get_compound(abbrev, lab_group=lab_group),
            "quantity":quantity,
            "unit":unit,
          })

      del dirty_data[reactant]

  # Apply the reactants to the clean data.
  for i, reactant in enumerate(reactants):
    clean_data["reactant_fk_"+str(i+1)] = reactant["reactant"]
    clean_data["quantity_"+str(i+1)] = reactant["quantity"]
    clean_data["unit_"+str(i+1)] = reactant["unit"]


  for key, val in clean_data.items():
    if not quick_validation(key, val):
      errors[key] = "Invalid value!"

  # Make sure the `ref` isn't duplicated.

  # Add any remaining field that don't require special processing.
  all_fields = get_model_field_names()
  for field in all_fields:
    if field not in clean_data and field in dirty_data:
      if quick_validation(field, dirty_data[field]):
        clean_data[field] = dirty_data[field]

  print "CLEANED___________________"
  print clean_data
  return (clean_data, errors)


def quick_validation(field, dirty_datum, model="Data"):
    if model=="Data":
        if dirty_datum=="?" and field in allow_unknown:
            return True
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



def clean_compound(compound):
  #Remove any non-printable characters.
  import string
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
  from DRP.models import get_lab_CG

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
