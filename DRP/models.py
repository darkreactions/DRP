from django.forms import *
from django.core import validators
from django.core.cache import cache

from django.contrib.auth.models import User, Group

from django.db import models
from django.db.models import Q

from validation import *
import random, string, datetime
from data_config import CONFIG
from chemspider_rdkit_extensions import *

#############  CACHE VALIDATION and ACCESS  ###########################
#Strip any spaces from the lab group title and/or the keys on cache access.
def set_cache(lab_group, key, value, duration=604800): #Default duration is 1 week.
 condensed_lab = lab_group.lab_title.replace(" ","")
 if key: condensed_key = key.replace(" ","") #Don't try to .replace None-types.
 else: condensed_key = None
 cache.set("{}|{}".format(condensed_lab, condensed_key), value, duration)

def get_cache(lab_group, key):
 condensed_lab = lab_group.lab_title.replace(" ","")
 condensed_key = key.replace(" ","") #Key must be a string.
 return cache.get("{}|{}".format(condensed_lab, condensed_key))

############### USER and LAB INTEGRATION #######################
ACCESS_CODE_MAX_LENGTH = 20 #Designates the max_length of access_codes

#Create a random alphanumeric code of specified length.
def get_random_code(length = ACCESS_CODE_MAX_LENGTH):
 return "".join(
   random.choice(
    string.letters + string.digits
   ) for i in range(length))


class Lab_Group(models.Model):
 lab_title = models.CharField(max_length=200)

 access_code = models.CharField(max_length=ACCESS_CODE_MAX_LENGTH,
  default=get_random_code)

 def __unicode__(self):
  return self.lab_title


############### USER CREATION #######################
class Lab_Member(models.Model):
 user = models.OneToOneField(User, unique=True) ###Allow lab member to switch?
 lab_group = models.ForeignKey(Lab_Group)

 def __unicode__(self):
  return self.user.username

############### COMPOUND GUIDE ########################
class CG_calculations(models.Model):
 json_data = models.TextField()
 compound = models.CharField(max_length=100, unique=True)
 smiles = models.CharField(max_length=200, unique=True)

 def __unicode__(self):
  return u"{} ({})".format(self.compound, self.smiles)

class CompoundEntry(models.Model):
 abbrev = models.CharField("Abbreviation", max_length=100) ###repr in admin 500 error
 compound = models.CharField("Compound", max_length=100)
 CAS_ID = models.CharField("CAS ID", max_length=13, blank=True)
 compound_type = models.CharField("Type", max_length=10)
 image_url = models.CharField("Image URL", max_length=100, blank=True)
 smiles = models.CharField("SMILES", max_length=100, blank=True)
 mw = models.CharField("Molecular Weight", max_length=20)

 lab_group = models.ForeignKey(Lab_Group, unique=False)
 #calculations = models.ForeignKey(CG_calculations)

 def __unicode__(self):
  if self.compound == self.abbrev:
   return u"{} (--> same) (LAB: {})".format(self.abbrev, self.lab_group.lab_title)
  return u"{} --> {} (LAB: {})".format(self.abbrev, self.compound, self.lab_group.lab_title)

TYPE_CHOICES = [[opt,opt] for opt in edit_choices["typeChoices"]]

def parse_CAS_ID(CAS):
  CAS = CAS.replace(" ", "-").replace("/", "-").replace("_", "-")

  #Check that the CAS ID has three hyphen-delineated parts.
  if len(CAS.split("-")) != 3:
   raise Exception("CAS ID requires three distinct parts.")
  #Check that only numbers are present.
  elif not CAS.replace("-","").isdigit():
   raise Exception("CAS ID may only have numeric characters.")

  return CAS


####errors[field] = "Field must be one of: {}".format(edit_choices[category])
def CG_validation(dirty_data, lab_group, editing_this=False):
 #Initialize the variables needed for the cleansing process.
 clean_data = {} #Keep track of cleaned fields
 errors = {}

 #Set default fields:
 clean_data["image_url"] = ""

 try:
  #If a CAS_ID is submitted, validate it.
  if dirty_data["CAS_ID"]:
   raw_CAS = dirty_data["CAS_ID"]
   clean_data["CAS_ID"] = parse_CAS_ID(raw_CAS) if raw_CAS else ""
  else:
   clean_data["CAS_ID"]=""
 except Exception as e:
  clean_data["CAS_ID"] = ""
  errors["CAS_ID"] = e

 other_fields = ["abbrev", "compound", "compound_type"]
 for field in other_fields:
  try:
   clean_data[field] = dirty_data[field]
  except:
   errors[field] = "This field cannot be blank."

 #Block duplicate abbrevs and compounds.
 if not editing_this:
  if not errors.get("compound") and CompoundEntry.objects.filter(compound=clean_data["compound"]).exists():
   errors["compound"] = "Compound already exists."
  if not errors.get("abbrev"):
   if CompoundEntry.objects.filter(abbrev=clean_data["abbrev"]).exists():
    errors["abbrev"] = "Abbreviation already used."

 #Make sure the compound exists in ChemSpider's database. ###Limits us to only ChemSpi?
 used_var = ""
 try:
  print clean_data["CAS_ID"]
  if clean_data["CAS_ID"]:
   used_var = "CAS_ID"
  elif not errors.get("compound"):
   used_var = "compound"
  else:
   return clean_data, errors

  query = chemspider_lookup({used_var:clean_data[used_var]})
  clean_data["image_url"], clean_data["smiles"], clean_data["mw"] = query
 except Exception as e:
  if used_var=="CAS_ID":
   errors["CAS_ID"] = "Could not validate CAS. Make sure it is correct."
  else:
   errors["compound"] = "Could not find this molecule. Try a different name."
 return clean_data, errors

######HERE


def update_all_compounds(lab_group=None):
 if lab_group:
  query = CompoundEntry.objects.filter(lab_group=lab_group)
  count = query.count()
  #Update all of the compounds before updating the reactions.
  print "Starting compound updates..."
  for i in query:
   update_compound(i, update_reactions=False)
  print "Starting reaction updates..."
  for i in query:
   update_reactions(i)
  print "Finished updating compounds for {}.".format(lab_group.lab_title)
 else:
  query = CompoundEntry.objects.all()
  print "Starting compound updates..."
  for i in query:
   update_compound(i, update_reactions=False)
  print "Starting reaction updates..."
  for i in query:
   update_reactions(i)
  print "Finished updating ALL compounds."

###REREAD TO PROVE USEFUL?
def collect_CG_entries(lab_group, overwrite=False):
 compound_guide = get_cache(lab_group, "COMPOUNDGUIDE")
 if not compound_guide or overwrite:
  compound_guide = CompoundEntry.objects.filter(lab_group=lab_group).order_by("compound")
  set_cache(lab_group, "COMPOUNDGUIDE", list(compound_guide))
 return compound_guide

def convert_QuerySet_to_list(query, model, with_headings=True):
 #Get the appropriate headings.
 all_fields = get_model_field_names(model=model, collect_ignored=True)
 fields_to_exclude = {"lab_group", "atoms"}
 headings = [field for field in all_fields if field not in fields_to_exclude]

 if with_headings:
  query_list = [list(headings)]
 else:
  query_list = []

 for entry in query:
  sub_list = []
  for field in headings:
   sub_list.append(getattr(entry, field))
  query_list.append(sub_list)

 return query_list

def collect_reactions_as_lists(lab_group, with_headings=True):
 try:
  if type(lab_group)==str:
   lab_group = Lab_Group.objects.filter(lab_title=lab_group)[0]
 except:
  raise Exception("Could not find lab group: \"{}\"".format(lab_group))

 query = Data.objects.filter(lab_group=lab_group)
 return convert_QuerySet_to_list(query, "Data", with_headings=with_headings)

def collect_CG_entries_as_lists(lab_group, with_headings=True):
 try:
  if type(lab_group)==str:
   lab_group = Lab_Group.objects.filter(lab_title=lab_group)[0]
 except:
  raise Exception("Could not find lab group: \"{}\"".format(lab_group))

 query = collect_CG_entries(lab_group)
 return convert_QuerySet_to_list(query, "CompoundEntry", with_headings=with_headings)

def collect_CG_name_pairs(lab_group, overwrite=False):
 pairs = get_cache(lab_group, "COMPOUNDGUIDE|NAMEPAIRS")
 if not pairs or overwrite:
  compound_guide = collect_CG_entries(lab_group)
  pairs = {entry.abbrev: entry.compound for entry in compound_guide}
  set_cache(lab_group, "COMPOUNDGUIDE|NAMEPAIRS", pairs)
 return pairs

def new_CG_entry(lab_group, **kwargs): ###Not re-read yet.
 try:
  new_entry = CompoundEntry()
  #Set the self-assigning fields:
  setattr(new_entry, "lab_group", lab_group)

  #Set the non-user field values.
  for (field, value) in kwargs.items(): #Assume data passed to the function is clean.
   setattr(new_entry, field, value)
  return new_entry
 except Exception as e:
  raise Exception("CompoundEntry construction failed!")

def get_good_rxns(lab_group=None, with_headings=True):
 #Convert string input into a lab_group if possible.
 try:
  if type(lab_group)==str:
    lab_group = Lab_Group.objects.filter(lab_title=lab_group)[0]

  #Collect ALL data or Lab-specific data.
  if lab_group:
   query = Data.objects.filter(lab_group=lab_group, is_valid=True)
  else:
   query = Data.objects.filter(is_valid=True)

  return convert_QuerySet_to_list(query, "Data", with_headings=with_headings)
 except:
  raise Exception("Cannot find lab_group: \"{}\"".format(lab_group))


############### DATA ENTRY ########################
calc_fields = ['XXXtitle', 'XXXinorg1', 'XXXinorg1mass', 'XXXinorg1moles', 'XXXinorg2', 'XXXinorg2mass',
   'XXXinorg2moles', 'XXXinorg3', 'XXXinorg3mass','XXXinorg3moles', 'XXXorg1', 'XXXorg1mass',
   'XXXorg1moles', 'XXXorg2', 'XXXorg2mass', 'XXXorg2moles', 'XXXoxlike1', 'XXXoxlike1mass',
   'XXXoxlike1moles', 'Temp_max', 'time', 'slowCool', 'pH',
   'leak', 'numberInorg', 'numberOrg', 'numberOxlike', 'numberComponents', 'inorgavgpolMax',
   'inorgrefractivityMax', 'inorgmaximalprojectionareaMax', 'inorgmaximalprojectionradiusMax',
   'inorgmaximalprojectionsizeMax', 'inorgminimalprojectionareaMax', 'inorgminimalprojectionradiusMax',
   'inorgminimalprojectionsizeMax', 'inorgavgpol_pHdependentMax',
   'inorgmolpolMax', 'inorgvanderwaalsMax', 'inorgASAMax',
   'inorgASA+Max', 'inorgASA-Max', 'inorgASA_HMax',
   'inorgASA_PMax', 'inorgpolarsurfaceareaMax', 'inorghbdamsaccMax', 'inorghbdamsdonMax',
   'inorgavgpolMin', 'inorgrefractivityMin', 'inorgmaximalprojectionareaMin',
   'inorgmaximalprojectionradiusMin', 'inorgmaximalprojectionsizeMin',
   'inorgminimalprojectionareaMin', 'inorgminimalprojectionradiusMin',
   'inorgminimalprojectionsizeMin', 'inorgavgpol_pHdependentMin',
   'inorgmolpolMin', 'inorgvanderwaalsMin', 'inorgASAMin',
   'inorgASA+Min', 'inorgASA-Min', 'inorgASA_HMin',
   'inorgASA_PMin', 'inorgpolarsurfaceareaMin', 'inorghbdamsaccMin', 'inorghbdamsdonMin',
   'inorgavgpolArithAvg', 'inorgrefractivityArithAvg', 'inorgmaximalprojectionareaArithAvg',
   'inorgmaximalprojectionradiusArithAvg', 'inorgmaximalprojectionsizeArithAvg',
   'inorgminimalprojectionareaArithAvg', 'inorgminimalprojectionradiusArithAvg',
   'inorgminimalprojectionsizeArithAvg',
   'inorgavgpol_pHdependentArithAvg', 'inorgmolpolArithAvg',
   'inorgvanderwaalsArithAvg', 'inorgASAArithAvg',
   'inorgASA+ArithAvg', 'inorgASA-ArithAvg',
   'inorgASA_HArithAvg', 'inorgASA_PArithAvg',
   'inorgpolarsurfaceareaArithAvg',
   'inorghbdamsaccArithAvg', 'inorghbdamsdonArithAvg',
   'inorgavgpolGeomAvg', 'inorgrefractivityGeomAvg',
   'inorgmaximalprojectionareaGeomAvg', 'inorgmaximalprojectionradiusGeomAvg',
   'inorgmaximalprojectionsizeGeomAvg', 'inorgminimalprojectionareaGeomAvg',
   'inorgminimalprojectionradiusGeomAvg', 'inorgminimalprojectionsizeGeomAvg',
   'inorgavgpol_pHdependentGeomAvg', 'inorgmolpolGeomAvg', 'inorgvanderwaalsGeomAvg',
   'inorgASAGeomAvg', 'inorgASA+GeomAvg', 'inorgASA-GeomAvg', 'inorgASA_HGeomAvg',
   'inorgASA_PGeomAvg', 'inorgpolarsurfaceareaGeomAvg', 'inorghbdamsaccGeomAvg', 'inorghbdamsdonGeomAvg',
   'orgavgpolMax', 'orgrefractivityMax',
   'orgmaximalprojectionareaMax', 'orgmaximalprojectionradiusMax', 'orgmaximalprojectionsizeMax',
   'orgminimalprojectionareaMax', 'orgminimalprojectionradiusMax', 'orgminimalprojectionsizeMax',
   'orgavgpol_pHdependentMax', 'orgmolpolMax',
   'orgvanderwaalsMax', 'orgASAMax', 'orgASA+Max', 'orgASA-Max', 'orgASA_HMax', 'orgASA_PMax',
   'orgpolarsurfaceareaMax', 'orghbdamsaccMax',
   'orghbdamsdonMax', 'orgavgpolMin', 'orgrefractivityMin',
   'orgmaximalprojectionareaMin', 'orgmaximalprojectionradiusMin', 'orgmaximalprojectionsizeMin',
   'orgminimalprojectionareaMin', 'orgminimalprojectionradiusMin',
   'orgminimalprojectionsizeMin', 'orgavgpol_pHdependentMin',
   'orgmolpolMin', 'orgvanderwaalsMin', 'orgASAMin',
   'orgASA+Min', 'orgASA-Min', 'orgASA_HMin', 'orgASA_PMin',
   'orgpolarsurfaceareaMin', 'orghbdamsaccMin',
   'orghbdamsdonMin', 'orgavgpolArithAvg', 'orgrefractivityArithAvg',
   'orgmaximalprojectionareaArithAvg', 'orgmaximalprojectionradiusArithAvg',
   'orgmaximalprojectionsizeArithAvg', 'orgminimalprojectionareaArithAvg',
   'orgminimalprojectionradiusArithAvg', 'orgminimalprojectionsizeArithAvg',
   'orgavgpol_pHdependentArithAvg', 'orgmolpolArithAvg', 'orgvanderwaalsArithAvg',
   'orgASAArithAvg', 'orgASA+ArithAvg', 'orgASA-ArithAvg', 'orgASA_HArithAvg', 'orgASA_PArithAvg',
   'orgpolarsurfaceareaArithAvg', 'orghbdamsaccArithAvg',
   'orghbdamsdonArithAvg', 'orgavgpolGeomAvg', 'orgrefractivityGeomAvg',
   'orgmaximalprojectionareaGeomAvg', 'orgmaximalprojectionradiusGeomAvg',
   'orgmaximalprojectionsizeGeomAvg', 'orgminimalprojectionareaGeomAvg',
   'orgminimalprojectionradiusGeomAvg', 'orgminimalprojectionsizeGeomAvg',
   'orgavgpol_pHdependentGeomAvg', 'orgmolpolGeomAvg', 'orgvanderwaalsGeomAvg',
   'orgASAGeomAvg', 'orgASA+GeomAvg', 'orgASA-GeomAvg',
   'orgASA_HGeomAvg', 'orgASA_PGeomAvg', 'orgpolarsurfaceareaGeomAvg', 'orghbdamsaccGeomAvg',
   'orghbdamsdonGeomAvg', 'oxlikeavgpolMax',
   'oxlikerefractivityMax', 'oxlikemaximalprojectionareaMax',
   'oxlikemaximalprojectionradiusMax', 'oxlikemaximalprojectionsizeMax', 'oxlikeminimalprojectionareaMax',
   'oxlikeminimalprojectionradiusMax', 'oxlikeminimalprojectionsizeMax',
   'oxlikeavgpol_pHdependentMax', 'oxlikemolpolMax',
   'oxlikevanderwaalsMax', 'oxlikeASAMax', 'oxlikeASA+Max',
   'oxlikeASA-Max', 'oxlikeASA_HMax', 'oxlikeASA_PMax',
   'oxlikepolarsurfaceareaMax', 'oxlikehbdamsaccMax',
   'oxlikehbdamsdonMax', 'oxlikeavgpolMin',
   'oxlikerefractivityMin', 'oxlikemaximalprojectionareaMin',
   'oxlikemaximalprojectionradiusMin', 'oxlikemaximalprojectionsizeMin',
   'oxlikeminimalprojectionareaMin', 'oxlikeminimalprojectionradiusMin',
   'oxlikeminimalprojectionsizeMin', 'oxlikeavgpol_pHdependentMin', 'oxlikemolpolMin',
   'oxlikevanderwaalsMin', 'oxlikeASAMin', 'oxlikeASA+Min',
   'oxlikeASA-Min', 'oxlikeASA_HMin', 'oxlikeASA_PMin',
   'oxlikepolarsurfaceareaMin', 'oxlikehbdamsaccMin',
   'oxlikehbdamsdonMin', 'oxlikeavgpolArithAvg',
   'oxlikerefractivityArithAvg', 'oxlikemaximalprojectionareaArithAvg',
   'oxlikemaximalprojectionradiusArithAvg', 'oxlikemaximalprojectionsizeArithAvg', 'oxlikeminimalprojectionareaArithAvg',
   'oxlikeminimalprojectionradiusArithAvg', 'oxlikeminimalprojectionsizeArithAvg', 'oxlikeavgpol_pHdependentArithAvg',
   'oxlikemolpolArithAvg', 'oxlikevanderwaalsArithAvg',
   'oxlikeASAArithAvg', 'oxlikeASA+ArithAvg',
   'oxlikeASA-ArithAvg', 'oxlikeASA_HArithAvg',
   'oxlikeASA_PArithAvg', 'oxlikepolarsurfaceareaArithAvg',
   'oxlikehbdamsaccArithAvg', 'oxlikehbdamsdonArithAvg',
   'oxlikeavgpolGeomAvg', 'oxlikerefractivityGeomAvg',
   'oxlikemaximalprojectionareaGeomAvg', 'oxlikemaximalprojectionradiusGeomAvg', 'oxlikemaximalprojectionsizeGeomAvg',
   'oxlikeminimalprojectionareaGeomAvg', 'oxlikeminimalprojectionradiusGeomAvg', 'oxlikeminimalprojectionsizeGeomAvg',
   'oxlikeavgpol_pHdependentGeomAvg', 'oxlikemolpolGeomAvg',
   'oxlikevanderwaalsGeomAvg', 'oxlikeASAGeomAvg',
   'oxlikeASA+GeomAvg', 'oxlikeASA-GeomAvg', 'oxlikeASA_HGeomAvg', 'oxlikeASA_PGeomAvg',
   'oxlikepolarsurfaceareaGeomAvg', 'oxlikehbdamsaccGeomAvg',
   'oxlikehbdamsdonGeomAvg', 'inorg-water-moleratio', 'inorgacc-waterdonratio', 'inorgdon-wateraccratio',
   'org-water-moleratio', 'orgacc-waterdonratio', 'orgdon-wateraccratio', 'inorg-org-moleratio',
   'inorgacc-orgdonratio', 'inorgdon-orgaccratio', 'notwater-water-moleratio', 'notwateracc-waterdonratio',
   'notwaterdon-wateraccratio', 'purity', 'outcome']

class DataCalc(models.Model):
 for calc_field in calc_fields:
  #Make sure field names don't contain operators.
  calc_field = calc_field.replace("+","PLUS").replace("-","MINUS")
  exec("{0} = models.CharField(\"{0}\", max_length=22)".format(calc_field))

 def __unicode__(self):
  return u"{}".format(self.XXXtitle);


#Create the form choices from the pre-defined ranges.
OUTCOME_CHOICES = [[opt,opt] for opt in edit_choices["outcomeChoices"]]
PURITY_CHOICES = [[opt,opt] for opt in edit_choices["purityChoices"]]
UNIT_CHOICES = [[opt,opt] for opt in edit_choices["unitChoices"]]
BOOL_CHOICES = [[opt,opt] for opt in edit_choices["boolChoices"]]

#Fields that are allowed to be stored as listy_strings.
list_fields = ["reactant", "quantity", "unit"]

#Many data are saved per lab group. Each data represents one submission.
class Data(models.Model):
 #List Fields
 for i in CONFIG.reactant_range():
  exec("reactant_{0} = models.CharField(\"Reactant {0}\", max_length=30)".format(i))
  exec("quantity_{0} = models.CharField(\"Quantity {0}\", max_length=10)".format(i))
  exec("unit_{0} = models.CharField(\"Unit {0}\", max_length=4)".format(i))

 ref = models.CharField("Reference", max_length=12)
 temp = models.CharField("Temperature", max_length=10)
 time = models.CharField("Time", max_length=10) ###
 pH = models.CharField("pH", max_length=5)

 #Yes/No/? Fields:
 slow_cool = models.CharField("Slow Cool", max_length=10)
 leak = models.CharField("Leak", max_length=10)
 outcome = models.CharField("Outcome", max_length=1)
 purity = models.CharField("Purity", max_length=1)

 notes = models.CharField("Notes", max_length=200, blank=True)

 #Self-assigning Fields:
 calculations = models.ForeignKey(DataCalc, unique=False, blank=True, null=True)
 calculated_pH = models.BooleanField()
 calculated_temp = models.BooleanField()
 calculated_time = models.BooleanField()

 atoms = models.CharField("Atoms", max_length=30, blank=True)

 user = models.ForeignKey(User, unique=False)
 lab_group = models.ForeignKey(Lab_Group, unique=False)
 creation_time = models.CharField("Created", max_length=26, null=True, blank=True)
 is_valid = models.BooleanField("Valid")
 #atoms = models.CharField("Atoms", max_length=40) ###Needs access to CG?

 #Categorizing Fields:
 duplicate_of = models.CharField("Duplicate", max_length=12, null=True, blank=True)
 recommended = models.CharField("Recommended", max_length=10)

 def __unicode__(self):
  return u"{} -- (LAB: {})".format(self.ref, self.lab_group.lab_title)

def validate_name(abbrev_to_check, lab_group):
 #Get the cached set of abbreviations.
 abbrevs = collect_CG_name_pairs(lab_group)
 return abbrev_to_check in abbrevs


#TODO:
def calculate_pH_from_CG(entry):
 pass

def calculate_pH_from_reaction(reaction_info):
 pass




#Add specified entries to a datum. Assume fields are now valid.
def new_Data_entry(user, **kwargs): ###Not re-read yet.
 try:
  new_entry = Data()
  #Set the self-assigning fields:
  setattr(new_entry, "creation_time", str(datetime.datetime.now()))
  setattr(new_entry, "user", user)
  setattr(new_entry, "lab_group", user.get_profile().lab_group)


  setattr(new_entry, "calculated_time", False)
  setattr(new_entry, "calculated_temp", False)

  #Set the non-user field values.
  for (field, value) in kwargs.items(): #Assume data passed to the function is clean.
   setattr(new_entry, field, value)

  #Calculate the pH if necessary.
  if not getattr(new_entry, "pH"):
   setattr(new_entry, "pH", calculate_pH_from_CG(new_entry))
   setattr(new_entry, "calculated_pH", True)
  else:
   setattr(new_entry, "calculated_pH", False)

  return new_entry
 except Exception as e:
  raise Exception("Data construction failed!")

def get_model_field_names(both=False, verbose=False, model="Data", unique_only=False, collect_ignored=False):
 clean_fields = []

 if model=="Data":
  if collect_ignored:
   fields_to_ignore = {u"id", "creation_time", "calculations"}
  else:
   fields_to_ignore = {u"id","user","lab_group", "atoms", "creation_time", "calculations", "calculated_temp", "calculated_time", "calculated_pH", "is_valid"}
  dirty_fields = [field for field in Data._meta.fields if field.name not in fields_to_ignore]
 elif model=="CompoundEntry":
  if collect_ignored:
   fields_to_ignore = {u"id", "image_url"}
  else:
   fields_to_ignore = {u"id","lab_group", "smiles", "mw"} ###Auto-populate?
  dirty_fields = [field for field in CompoundEntry._meta.fields if field.name not in fields_to_ignore]
 else:
  raise Exception("Unknown model specified.")

 #Ignore any field that is in fields_to_ignore.
 for field in dirty_fields:
  #Return the non list-fields:
  if unique_only and field.name[-1].isdigit(): continue

  #Return either the verbose names or the non-verbose names.
  if both:
   clean_fields += [{"verbose":field.verbose_name, "raw":field.name}] ###Make verbose names pretty
  elif verbose:
   clean_fields += [field.verbose_name] ###Make verbose names pretty
  else:
   clean_fields += [field.name]
 return clean_fields

def get_ref_set(lab_group):
 ref_set = get_cache(lab_group, "REFS")
 if not ref_set:
  ref_set = set(Data.objects.values_list('ref', flat=True))
  set_cache(lab_group, "REFS", ref_set)
 return ref_set

def revalidate_data(data, lab_group, batch=False):
 #Collect the data to validate
 dirty_data = {field:getattr(data, field) for field in get_model_field_names()}
 #Validate and collect any errors
 (clean_data, errors) = full_validation(dirty_data, lab_group, revalidating=True)

 is_valid = False if errors else True

 if is_valid:
  print "VALIDATED A BAD DATUM!"

 setattr(data, "is_valid", is_valid)
 data.save()

 return errors

 #Does not auto-clear the cache --only modifies the data entry.

def full_validation(dirty_data, lab_group, revalidating=False):
 parsed_data = {} #Data that needs to be checked.
 clean_data = {} #Keep track of cleaned fields
 errors = {}

 fields = get_model_field_names()

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
     clean_data[field] = "" #If nothing was entered, store nothing ###Used to be "?" -- why?
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

    #Check to make sure no references are repeated.
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
