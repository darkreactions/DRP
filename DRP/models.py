from django.forms import *
from django.core import validators
from django.core.cache import cache
from django.contrib.auth.models import User, Group
from django.db import models
from django.db.models import Q

from data_config import CONFIG
from validation import *
from uuid import uuid4
from CGCalculator import CGCalculator
from collections import defaultdict
from subprocess import Popen
from DRP.settings import LOG_DIR, BASE_DIR

import json, random, string, datetime, operator
import rdkit.Chem as Chem
import chemspipy


#Basic Retrieval Functions Necessary in Models:
#Get the data that belongs to a Lab_Group
def get_lab_Data(lab_group):
 return Data.objects.filter(lab_group=lab_group).order_by("creation_time_dt")

def get_ref_set(lab_group, reset_cache=True):
 ref_set = get_cache(lab_group, "DATAREFS")
 if not ref_set or reset_cache:
  ref_set = set(get_lab_Data(lab_group).values_list('ref', flat=True))
  set_cache(lab_group, "DATAREFS", ref_set)
 return ref_set

def get_Lab_Group(query):
 try:
  if type(query==Lab_Group):
   return query
  return Lab_Group.objects.filter(lab_title=raw_string).get()
 except:
  raise Exception("Could not find Lab_Group with lab_title: {}".format(raw_string))

def get_lab_CG(lab_group):
 return CompoundEntry.objects.filter(lab_group=lab_group).order_by("compound")


# # # # # # # # # # # # # # # # # # #
  # # # # # # # # RDKIT and ChemSpider Functions # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # #
def get_first_chemspider_entry(search_fields):
 for i in search_fields:
  try:
   query = chemspipy.find_one(i)
   return query
  except Exception as e:
   pass
 return None

def chemspider_lookup(val):
 chemspi_query = ""
 #Accept either a CompoundEntry object or a dict with the valid fields
 search_fields = ["CAS_ID", "compound"]
 if type(val)==dict:
  query_criteria = [val.get(i) for i in search_fields if val.get(i)]
 elif type(val)==CompoundEntry:
  query_criteria = [getattr(val, i) for i in search_fields]
 else:
  query_criteria = [val]

 try:
  query = get_first_chemspider_entry(search_fields)
  assert query
  return query
 except:
  raise Exception("Could not find compound on ChemSpider!")

def get_atoms_from_compound(CG_entry = None):
 return get_atoms_from_smiles(CG_entry.smiles)

def get_atoms_from_smiles(smiles, show_hydrogen=False):
 if not smiles:
  raise Exception("SMILES cannot be None!")

 mols = Chem.MolFromSmiles(str(smiles),sanitize=False)
 if mols == None:
  return []
 #TODO: Incorporate hydrogens into model.
 #if show_hydrogen:
 # try:
 #  mols = Chem.AddHs(mols) ###PRECONDITION?
 # except:
 #  pass
 atoms = mols.GetAtoms()
 return [atom.GetSymbol() for atom in atoms]

def get_atom_count_from_smiles(smiles):
 atom_list = get_atoms_from_smiles(smiles)

 #Count the number of occurances of each atom.
 atom_count = defaultdict(int)
 for atom in atom_list:
  atom_count[atom] += 1

 return dict(atom_count) #Return a normal dictionary, not a defaultdict.

def collect_CGs_by_abbrevs(lab_group, abbrev_list):
 CG_list = []
 for i in abbrev_list:
  query = get_lab_CG(lab_group).filter(abbrev=i)
  if query.exists():
   CG_list.append(query[0])
 return CG_list

def get_smiles_from_CG_list(CG_list, allow_custom=True):
 smiles_list = [i.smiles for i in CG_list if (not i.custom or allow_custom)]
 return smiles_list

def condense_smiles_list_to_atoms(smiles_list):
 atoms_list = []
 for i in smiles_list:
  if i:
   atoms_list += get_atoms_from_smiles(i)
 return set(atoms_list)

def get_abbrevs_from_reaction(reaction):
 abbrevs_list = [getattr(reaction, "reactant_{}".format(i)) for i in CONFIG.reactant_range() if getattr(reaction, "reactant_{}".format(i))]
 return abbrevs_list

def get_atom_set_from_abbrevs(lab_group, abbrev_list):
 return condense_smiles_list_to_atoms(
  get_smiles_from_CG_list(
   collect_CGs_by_abbrevs(lab_group, abbrev_list),
   allow_custom=False
   ))

def get_atom_set_from_reaction(reaction):
 return get_atom_set_from_abbrevs(reaction.lab_group, get_abbrevs_from_reaction(reaction))

def update_compound_and_reactions(lab_group, entry):
 try:
  update_compound(entry)
  update_reactions_with_compound(lab_group, entry)
 except Exception as e:
  print e
  raise Exception("Compound_and_reactions update failed!")

def refresh_compound_guide(lab_group = None, verbose = False):
 #Either refresh all of the data or the data for a specific lab group.
 if lab_group:
  query = get_lab_CG(lab_group)
 else:
  query = CompoundEntry.objects.all()

 #Actually perform the refresh/update.
 i = 0
 for compound in query:
  try:
   if verbose:
    i += 1
    if i%10==0: print "... {}.".format(i)
   update_compound(compound)
  except Exception as e:
   print "Could not update: {}".format(compound)

#Update the compound by reloading the ChemSpider search data.
def update_compound(entry):
  try:
    if not entry.custom: #Only update compounds that are not custom.
      #Get the most up-to-date ChemSpider info for a given CAS/compound.
      query = get_first_chemspider_entry([entry.CAS_ID, entry.compound])
      if query:
        #Update the entry.
        entry.image_url, entry.smiles, entry.mw = query.imageurl, query.smiles, query.molecularweight
        perform_calcs = True

      else:
        perform_calcs = False
        print "Found legacy entry that should be custom: {}".format(entry.compound)
    else:
        perform_calcs = False
        entry.calculations = None
        entry.image_url, entry.smiles, entry.mw = "","",""
    entry.save()

    if perform_calcs:
      #Start a new compound-calc worker process to determine compound properties.
      err_log = open(LOG_DIR+"/compound_calculations/error.log","a")
      act_log = open(LOG_DIR+"/compound_calculations/process.log","a")
      worker_script = BASE_DIR+"/DRP/compound_calculations/calculate_compound_properties.py"
      command = "python {} {}".format(worker_script, entry.id)
      #Log to the files above and make the worker independent of the parent process.
      Popen(command.split(), stdout=act_log, stderr=err_log, close_fds=True)

  except Exception as e:
    print e
    raise Exception("Compound update failed!")

def update_reaction(reaction, lab_group):
 #Store the atoms as a string -- not a set.
 reaction.atoms = "".join(get_atom_set_from_reaction(reaction))
 #Revalidate and save the datum.
 revalidate_datum(reaction, lab_group)

def update_reactions_with_compound(lab_group, compound):
 #Update the individual "atom" records on each reaction.
 lab_data = get_lab_Data(lab_group)
 changed_reactions = get_Data_with_abbrev(lab_data, compound)
 for reaction in changed_reactions:
  try:
   update_reaction(reaction, lab_group)
  except Exception as e:
   print "update_reactions_with_compound failed: {}".format(e)

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
 lab_title = models.CharField(max_length=200, unique=True, error_messages={'unique':"This name is already taken."})
 lab_address = models.CharField(max_length=200)
 lab_email = models.CharField(max_length=254) #Maximum length of email address
 access_code = models.CharField(max_length=ACCESS_CODE_MAX_LENGTH,
  default=get_random_code)

 def __unicode__(self):
  return self.lab_title

############### USER CREATION #######################
class Lab_Member(models.Model):
 user = models.OneToOneField(User, related_name="profile", unique=True)
 license_agreement_date_dt = models.DateTimeField("License Agreement Date", null=True, blank=True)
 lab_group = models.ForeignKey(Lab_Group)

 def update_license(self):
  self.license_agreement_date_dt = datetime.datetime.now()
  self.save()

 def __unicode__(self):
  return self.user.username

############### DATA ENTRY ########################
class DataCalc(models.Model):
  contents = models.TextField(default="[]")

  def __init__(self, jsonContent):
    self.contents = json.dumps(jsonContent)

  def __unicode__(self):
    return u"{}".format(self.contents);

  #Convert the stringy contents to an actual array/JSON object.
  def make_json(self):
    return json.loads(self.contents)


#Many data are saved per lab group. Each data represents one submission.
class Data(models.Model):
 ref = models.CharField("Reference", max_length=12)

 #List Fields
 for i in CONFIG.reactant_range():
  exec("reactant_{0} = models.CharField(\"Reactant {0}\", max_length=30)".format(i))
  exec("quantity_{0} = models.CharField(\"Quantity {0}\", max_length=10)".format(i))
  exec("unit_{0} = models.CharField(\"Unit {0}\", max_length=4)".format(i))

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
 calculated_pH = models.BooleanField(default=False)
 calculated_temp = models.BooleanField(default=False)
 calculated_time = models.BooleanField(default=False)

 atoms = models.CharField("Atoms", max_length=30, blank=True)

 user = models.ForeignKey(User, unique=False)
 lab_group = models.ForeignKey(Lab_Group, unique=False)
 creation_time_dt = models.DateTimeField("Created", null=True, blank=True)
 is_valid = models.BooleanField("Valid", default=False)

 #Categorizing Fields:
 public = models.BooleanField("Public", default=False)
 duplicate_of = models.CharField("Duplicate", max_length=12, null=True, blank=True)
 recommended = models.CharField("Recommended", max_length=10)

 def __unicode__(self):
  return u"{} -- (LAB: {})".format(self.ref, self.lab_group.lab_title)

############### RECOMMENDATIONS ########################
class ModelStats(models.Model):
  # model false-positive on test set
  false_positive_rate = models.FloatField()

  # recommendation quality
  actual_success_rate = models.FloatField()

  # evaluation of similarity metric
  estimated_success_rate = models.FloatField()

  # model performance
  performance = models.FloatField()
  datetime = models.DateTimeField()

  #Model Descriptors
  title = models.CharField("Title", max_length=100, default="")
  description = models.TextField(default="")

  def __unicode__(self):
    return "Performance:{} ({})".format(self.performance, self.datetime)

class Model_Version(models.Model):
 model_type = models.CharField("Type", max_length=20)
 date_dt = models.DateTimeField("Date", null=True, blank=True)
 notes = models.CharField("Notes", max_length=200, blank=True)
 lab_group = models.ForeignKey(Lab_Group)

class Recommendation(models.Model):
 #Reactant Fields
 for i in CONFIG.reactant_range():
  exec("reactant_{0} = models.CharField(\"Reactant {0}\", max_length=30)".format(i))
  exec("quantity_{0} = models.CharField(\"Quantity {0}\", max_length=10)".format(i))
  exec("unit_{0} = models.CharField(\"Unit {0}\", max_length=4)".format(i))

 score = models.FloatField("Score")
 temp = models.CharField("Temperature", max_length=10)
 time = models.CharField("Time", max_length=10) ###
 pH = models.CharField("pH", max_length=5)

 #Yes/No/? Fields:
 slow_cool = models.CharField("Slow Cool", max_length=10)
 leak = models.CharField("Leak", max_length=10)
 outcome = models.CharField("Outcome", max_length=1)
 purity = models.CharField("Purity", max_length=1)

 #Self-assigning Fields:
 atoms = models.CharField("Atoms", max_length=30, blank=True)
 lab_group = models.ForeignKey(Lab_Group, unique=False)
 model_version = models.ForeignKey(Model_Version, unique=False)
 user = models.ForeignKey(User, unique=False, null=True, blank=True, default=None, related_name="last_user")
 assigned_user = models.ForeignKey(User, unique=False, null=True, blank=True, default=None, related_name="assigned_user")
 seed = models.ForeignKey(Data, unique=False, null=True, blank=True, default=None)
 date_dt = models.DateTimeField("Created", null=True, blank=True)
 complete = models.BooleanField("Complete", default=False)
 seeded = models.BooleanField("From Seed", default=False)

 #Fields for user feedback.
 saved = models.BooleanField("Saved", default=False)
 nonsense = models.BooleanField("Nonsense", default=False)
 hidden = models.BooleanField("Hidden", default=False)
 notes = models.CharField("Notes", max_length=200, blank=True)

 def __unicode__(self):
  return u"REC: {} -- (LAB: {} -- Saved: {})".format(self.score, self.lab_group.lab_title, self.saved)

class RankedReactionList(models.Model):
 original_list = models.TextField()
 seed = models.TextField()
 ranked_list = models.TextField()
 ranker = models.ForeignKey(User, unique=False, default=None, null=True)

 def __unicode__(self):
  return u"RANKED_LIST: Seed: {} -- (Ranker: {})".format(self.seed, self.ranker)

 def get_original_list(self):
  return json.loads(self.original_list)

 def get_ranked_list(self):
  return json.loads(self.ranked_list)

 def get_seed(self):
  return json.loads(self.seed)

 #Returns a shuffled list and the index of the seed.
 def get_shuffled_list(self):
  shuffled_list = self.get_original_list()
  random.shuffle(shuffled_list)
  return (shuffled_list)

 def store_ranked_list(self, ranked_list, ranker):
  self.ranker = ranker
  self.ranked_list = json.dumps(ranked_list)
  self.save()

def get_unranked_reactions(seed=None):
 unranked = RankedReactionList.objects.filter(ranker=None)
 #If a seed is specified, apply it to the filter.
 if seed:
  #Convert any lists to strings to allow filtering.
  seed = json.dumps(seed) if type(seed)==list else seed
  unranked = unranked.filter(seed=seed)
 return unranked

def get_random_unranked_reaction_or_none(seed=None):
 unranked_rxns = get_unranked_reactions(seed=seed)
 if unranked_rxns.exists():
  random_index = random.randrange(unranked_rxns.count())
  return unranked_rxns[random_index]
 return None

############### COMPOUND GUIDE ########################
class CG_calculations(models.Model):
 json_data = models.TextField()
 compound = models.CharField(max_length=200)
 smiles = models.CharField(max_length=255, unique=True)
 json = models.TextField(null=True, default="{}")

 def __unicode__(self):
  return u"{} ({})".format(self.compound, self.smiles)

def create_CG_calcs_if_needed(compound, smiles, compound_type):
    #Variable Setup
    jchem_path =  CONFIG.jchem_path
    sdf_path = "tmp"

    compound, smiles, compound_type = str(compound), str(smiles), str(compound_type)

    #Only Organics that have smiles may have calculations.
    if compound_type not in {"Org", "Inorg"} or not smiles:
        return
    #Either return an old CG_calculation or a new one.
    print compound_type
    try:
        cgc = CG_calculations.objects.filter(smiles=smiles)[0]
    except Exception as e:
        #Calculate properties for the CGEntry
        sdf_filename = str(uuid4()) + filter(str.isalnum, compound)
        #TODO: Speed this up? This is dreadfully slow.
        props = CGCalculator(compound, sdf_filename, smiles, compound_type, jchem_path, sdf_path).get_properties()
        props = json.dumps(props)
        #Store the actual CG_calculation in the database.
        cgc = CG_calculations(json_data=props, compound=compound, smiles=smiles)
        cgc.save()

	#Set the calculations field in each CompoundEntry.
	CompoundEntry.objects.filter(smiles=smiles).update(calculations=cgc)
    return cgc

def perform_CG_calculations(only_missing=True, lab_group=None, reattempt_failed = False, verbose=False):
 #Variable Setup
 success = 0
 i = 0

 cg = CompoundEntry.objects.all()
 if only_missing:
  cg = cg.filter(calculations=None)
 if lab_group:
  cg = cg.filter(lab_group=lab_group)
 if not reattempt_failed:
  cg = cg.filter(calculations_failed=False)

 for entry in cg:
  try:
   if verbose:
    i+=1
    if i%5==0: print "... {}.".format(i)

   try:
    entry.compound = clean_compound(entry.compound)
    calc = create_CG_calcs_if_needed(entry.compound, entry.smiles, entry.compound_type)
    entry.calculations = calc
   except:
    entry.calculations_failed = True

   entry.save
   success += 1
  except Exception as e:
   print "CG_calculation construction failed: {}".format(entry.compound)
   print "ERROR: {}".format(e)

 print "CG_calculations complete! ({} of {} entries changed)".format(success, cg.count())

#Remove any non-printable characters.
def clean_compound(compound):
 return filter(lambda x: x in string.printable, compound)

class CompoundEntry(models.Model):
 abbrev = models.CharField("Abbreviation", max_length=100)
 compound = models.CharField("Compound", max_length=100)
 CAS_ID = models.CharField("CAS ID", max_length=13, blank=True, default="")
 compound_type = models.CharField("Type", max_length=10)
 image_url = models.CharField("Image URL", max_length=100, blank=True, default="")
 smiles = models.CharField("SMILES", max_length=255, blank=True, default="")
 mw = models.CharField("Molecular Weight", max_length=20, default="")
 custom = models.BooleanField("Custom", default=False)

 lab_group = models.ForeignKey(Lab_Group, unique=False)
 calculations = models.ForeignKey(CG_calculations, unique=False, null=True, default=None)
 calculations_failed = models.BooleanField(default=False)

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

def validate_CG(dirty_data, lab_group, editing_this=False):
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
   query = get_first_chemspider_entry(search_fields)
   clean_data["image_url"], clean_data["smiles"], clean_data["mw"] = query.imageurl, query.smiles, query.molecularweight
  except:
   if search_fields[0]:
    errors["CAS_ID"] = "Could not find a molecule with this CAS ID."
   else:
    errors["compound"] = "Could not find this compound."

 #Prevent duplicate abbrevs.

 if not errors.get("abbrev"):
  print errors
  clean_data["abbrev"] = dirty_data["abbrev"]
  if not editing_this:
   if get_lab_CG(lab_group).filter(abbrev=clean_data["abbrev"]).exists():
    errors["abbrev"] = "Abbreviation already used."
 return clean_data, errors


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

def convert_Data_to_list(dat, headings=None):
	if not headings:
		all_fields = get_model_field_names(model="Data", collect_ignored = True)
		fields_to_exclude = {"lab_group", "atoms"}
		headings = [field for field in all_fields if field not in fields_to_exclude]

	r_heads = ["reactant_" + str(i) for i in range(1,6)]

	results = []
	for field in headings:
		val = getattr(dat, field)
		if field in r_heads and val != '':
			new_val = CompoundEntry.objects.filter(abbrev=val).first()
			if new_val is None:
				new_val = CompoundEntry.objects.filter(compound=val).first()
			new_val = new_val.compound
			if new_val:
				val = new_val
		results.append(val)

	return results

def collect_reactions_as_lists(lab_group, with_headings=True):
 lab_group = get_Lab_Group(lab_group)

 query = get_lab_Data(lab_group)
 return convert_QuerySet_to_list(query, "Data", with_headings=with_headings)

def get_lab_CG_as_lists(lab_group, with_headings=True):
 lab_group = get_Lab_Group(lab_group)

 query = get_lab_CG(lab_group)
 return convert_QuerySet_to_list(query, "CompoundEntry", with_headings=with_headings)

def collect_CG_name_pairs(lab_group, overwrite=False):
 pairs = get_cache(lab_group, "COMPOUNDGUIDE|NAMEPAIRS")
 if not pairs or overwrite:
  compound_guide = get_lab_CG(lab_group)
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

#Filter the Data by a specific abbrev.
def get_Data_with_abbrev(lab_data, abbrev):
 if type(abbrev)==CompoundEntry:
  abbrev = abbrev.abbrev #Meta...

 Q_list = [Q(("reactant_{}".format(i),abbrev)) for i in CONFIG.reactant_range()]
 return lab_data.filter(reduce(operator.or_, Q_list))

#Collect a list of all valid data either globally or for a specific lab.
def get_good_rxns(lab_group=None, with_headings=True):
 if lab_group:
  lab_group = get_Lab_Group(lab_group)
  query = get_lab_Data(lab_group).filter(is_valid=True)
 else:
   query = Data.objects.filter(is_valid=True)

 return convert_QuerySet_to_list(query, "Data", with_headings=with_headings)



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
  setattr(new_entry, "creation_time_dt", datetime.datetime.now())
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

def get_model_field_names(both=False, verbose=False, model="Data", unique_only=False, collect_ignored=False, for_upload=False):
 clean_fields = []

 if model=="Data":
  if collect_ignored:
   fields_to_ignore = {u"id", "creation_time_dt", "calculations"}
  else:
   fields_to_ignore = {u"id","user","lab_group", "atoms", "creation_time_dt",
                       "calculations", "calculated_temp", "calculated_time",
                       "calculated_pH", "is_valid", "public"}
  dirty_fields = [field for field in Data._meta.fields if field.name not in fields_to_ignore]
 elif model=="Recommendation":
  if collect_ignored:
   fields_to_ignore = {u"id", "creation_time_dt"}
  else:
   fields_to_ignore = {u"id","user", "assigned_user", "lab_group", "saved",
                       "model_version", "atoms", "creation_time_dt", "nonsense",
                       "complete", "score", "date_dt", "hidden", "seed", "seeded"}
  dirty_fields = [field for field in Recommendation._meta.fields if field.name not in fields_to_ignore]
 elif model=="CompoundEntry":
  if collect_ignored:
   fields_to_ignore = {u"id", "image_url", "custom", "calculations"}
  else:
   fields_to_ignore = {u"id","lab_group", "smiles", "mw", "custom",
                       "calculations", "calculations_failed"}
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

######################  Developer Functions  ###########################

def update_all_reactions(lab_group):
 data = get_lab_Data(lab_group)
 for entry in data:
  try:
   update_reaction(reaction, lab_group)
  except:
   print "--could not update reaction: {}".format(entry)
 print "Finished data validation."
