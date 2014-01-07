from django.forms import *
from django.core import validators
from django.core.cache import cache

from django.contrib.auth.models import User, Group

from django.db import models
from django.db.models import Q

from validation import *
import random, string, datetime
from data_config import CONFIG
#from chemspider_rdkit_extensions import *
import rdkit.Chem as Chem
import chemspipy

# # # # # # # # # # # # # # # # # # #
  # # # # # # # # RDKIT and ChemSpider Functions # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # #
def chemspider_lookup(cg_entry):
 chemspi_query = ""
 try:
  #Accept either a CompoundEntry object or a dict with the valid fields
  search_fields = ["CAS_ID", "compound"]
  if type(cg_entry)==dict:
   query_criteria = [cg_entry.get(i) for i in search_fields if cg_entry.get(i)]
  elif type(cg_entry)==CompoundEntry:
   query_criteria = [getattr(cg_entry, i) for i in search_fields]
  elif type(cg_entry)==str:
   query_criteria = [cg_entry]

  #Works for both dicts and CG_entries
  assert(query_criteria)
 except:
  raise Exception("chemspider_lookup accepts strings, dicts, or CompoundEntries")

 for i in query_criteria:
  if i: #ChemSpider doesn't like empty requests.
   chemspi_query = chemspipy.find_one(i)
  if chemspi_query: break

 if chemspi_query:
  return chemspi_query.imageurl, chemspi_query.smiles, chemspi_query.molecularweight
 else:
  raise Exception("Could not find compound on ChemSpider!")

def get_atoms_from_compound(CG_entry = None):
 return get_atoms_from_smiles(CG_entry.smiles)

def get_atoms_from_smiles(smiles, show_hydrogen=False):
 if not smiles:
  print "No smiles found!"
  return []

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

def collect_CGs_by_abbrevs(lab_group, abbrev_list):
 CG_list = []
 for i in abbrev_list:
  query = CompoundEntry.objects.filter(lab_group=lab_group, abbrev=i)
  if query.exists(): #Don't append index an empty query; Django gets annoyed.
   CG_list.append(query[0])
 return CG_list

def get_smiles_from_CG_list(CG_list):
 smiles_list = [i.smiles for i in CG_list]
 return smiles_list

def condense_smiles_list_to_atoms(smiles_list, show_hydrogen = False):
 atoms_list = []
 for i in smiles_list:
  atoms_list += get_atoms_from_smiles(i, show_hydrogen)
 return set(atoms_list)

def get_abbrevs_from_reaction(reaction):
 abbrevs_list = [getattr(reaction, "reactant_{}".format(i)) for i in CONFIG.reactant_range() if getattr(reaction, "reactant_{}".format(i))]
 return abbrevs_list

def get_atom_set_from_abbrevs(lab_group, abbrev_list):
 return condense_smiles_list_to_atoms(
  get_smiles_from_CG_list(
   collect_CGs_by_abbrevs(lab_group, abbrev_list)
   ), show_hydrogen=True
  )

def get_atom_set_from_reaction(reaction):
 return get_atom_set_from_abbrevs(reaction.lab_group, get_abbrevs_from_reaction(reaction))

def get_reactions_with_compound(compound):
 Q_string = ""
 abbrev = compound.abbrev
 for i in CONFIG.reactant_range():
  Q_string += "Q(reactant_{}=\"{}\")|".format(i, abbrev)

 if Q_string:
  Q_string = Q_string[:-1] #Remove the last trailing "|" if applicable.

 return eval("Data.objects.filter(lab_group=compound.lab_group).filter({})".format(Q_string))

def update_compound(compound, update_data=True):
 try:
  #Update the CG entry itself. Make sure "inorg" types don't query ChemSpider.
  try:
   compound.image_url, compound.smiles, compound.mw = chemspider_lookup(compound)
   print compound.compound_type, compound.abbrev
  except:
   if compound.compound_type=="Inorg":
    compound.image_url, compound.smiles, compound.mw = "","",""
   else:
    raise Exception("Could not find via ChemSpider!")
  compound.save()

  #Update the individual "atom" records on each reaction.
  if update_data:
   update_reactions(compound)
 except Exception as e:
  print "Could not update {}\n\t{}".format(compound, e)

def update_reactions(compound):
 #Update the individual "atom" records on each reaction.
 changed_reactions = get_reactions_with_compound(compound)
 for reaction in changed_reactions:
  #Store the atoms as a string -- not a set.
  reaction.atoms = "".join(get_atom_set_from_reaction(reaction))
  reaction.save()

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

def get_Lab_Group(raw_string):
 try:
  Lab_Group.objects.filter(lab_title=raw_string)
 except:
  raise Exception("Could not find Lab_Group with lab_title: {}".format(raw_string))

############### USER CREATION #######################
class Lab_Member(models.Model):
 user = models.OneToOneField(User, unique=True)
 license_agreement_date = models.CharField("License Agreement Date", max_length=26, blank=True) ###TODO: Explore why this isn't a datetime field. 
 lab_group = models.ForeignKey(Lab_Group)

 def is_licensed(self):
  return self.license_agreement_date and (self.license_agreement_date > CONFIG.license_agreement_date)

 def __unicode__(self):
  return self.user.username

############### RECOMMENDATIONS ########################
class Model_Version(models.Model):
 model_type = models.CharField("Type", max_length=20)
 date = models.CharField("Date", max_length=26) ###TODO: Explore why this isn't a datetime field. 
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
 date = models.CharField("Created", max_length=26, null=True, blank=True) ###TODO: Explore why this isn't a datetime field. 

 #Fields for user feedback.
 saved = models.BooleanField("Saved", default=False)
 nonsense = models.BooleanField("Nonsense", default=False)
 notes = models.CharField("Notes", max_length=200, blank=True)

 def __unicode__(self):
  return u"REC: score -- (LAB: {} -- Saved: {})".format(self.score, self.lab_group.lab_title, self.saved)

#Organize the reactants into sublists: [[reactant, quantity, unit], ...]
def partition_reactant_fields(lst):
 num_fields = CONFIG.fields_per_reactants
 total_fields = num_fields * CONFIG.num_reactants
 #Sort the <total_fields> elements of lst into sublists of <num_fields> elements each
 reaction = []
 reactant = []
 for (datum, i) in zip(lst, xrange(total_fields)):
  if i%num_fields==0 and i!=0:
   reaction.append(reactant)
   reactant = [datum] 
  else:
   reactant.append(datum)
 reaction.append(reactant)
 return reaction
 
#Apply the headings to the variables in the list.
###TODO: Notes for Paul: Remove the "ref"
def get_vars_from_list(lst):
 verbose_headers = get_model_field_names(verbose=True, model="Data")[16:]
 info = [[i, j] for (i,j) in zip(verbose_headers, lst)]
 return info 

#Creates the Recommendation entry, but does not store it in database.
def field_list_to_Recommendation(lab_group, lst, in_bulk=False):
 try:
  new_rec = Recommendation()
  #Set the self-assigning fields:
  setattr(new_rec, "lab_group", lab_group)
  setattr(new_rec, "score", float(lst[0]))

  #Set the non-user field values.
  fields = get_model_field_names(model="Recommendation")
  for (field, value) in zip(fields, lst[2][1:]): #Ignore the reference field.
   setattr(new_rec, field, value)
  return new_rec
   
  new_rec.atoms = "".join(get_atom_set_from_abbrevs(lab_group, get_abbrevs_from_reaction(new_rec)))
 
  if not in_bulk:
    new_rec.date = str(datetime.datetime.now())

  return new_rec
 except Exception as e:
  raise Exception("Recommendation construction failed: {}".format(e))

def store_new_Recommendation_list(lab_group, list_of_recommendations, version_notes = ""):
 if type(lab_group) != Lab_Group:
  lab_group = get_Lab_Group(lab_group)
 
 call_time = str(datetime.datetime.now())
 
 #Store the information for this "Version" of the Recommendation model.
 new_version = Model_Version()
 new_version.model_type = "Recommendation"
 new_version.date = call_time
 new_version.notes = version_notes
 new_version.lab_group = lab_group
 new_version.save()

 #Store the actual Recommendation entries.
 num_success = 0
 count = 0
 for i in list_of_recommendations:
  count += 1
  try:
   new_rec = field_list_to_Recommendation(lab_group, i, in_bulk=True)
   new_rec.date = call_time
   new_rec.model_version = new_version
   new_rec.save() #Store this recommendation in the database
   num_success += 1
  except Exception as e:
    print "Recommendation {} could not be constructed: {}".format(count, e)
 
 print "Finished creating and storing {} of {} items!.".format(num_success, count)

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
 if type(lab_group) != Lab_Group:
  lab_group = get_Lab_Group(lab_group)

 query = Data.objects.filter(lab_group=lab_group)
 return convert_QuerySet_to_list(query, "Data", with_headings=with_headings)

def collect_CG_entries_as_lists(lab_group, with_headings=True):
 if type(lab_group) != Lab_Group:
  lab_group = get_Lab_Group(lab_group)

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
 if lab_group:
  if type(lab_group) != Lab_Group:
   lab_group = get_Lab_Group(lab_group)
  query = Data.objects.filter(lab_group=lab_group, is_valid=True)
 else:
   query = Data.objects.filter(is_valid=True)

 return convert_QuerySet_to_list(query, "Data", with_headings=with_headings)


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
 calculated_pH = models.BooleanField(default=False)
 calculated_temp = models.BooleanField(default=False)
 calculated_time = models.BooleanField(default=False)

 atoms = models.CharField("Atoms", max_length=30, blank=True)

 user = models.ForeignKey(User, unique=False)
 lab_group = models.ForeignKey(Lab_Group, unique=False)
 creation_time = models.CharField("Created", max_length=26, null=True, blank=True)
 is_valid = models.BooleanField("Valid", default=False)

 #Categorizing Fields:
 public = models.BooleanField("Public", default=False)
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
   fields_to_ignore = {u"id","user","lab_group", "atoms", "creation_time", "calculations", "calculated_temp", "calculated_time", "calculated_pH", "is_valid", "public"}
  dirty_fields = [field for field in Data._meta.fields if field.name not in fields_to_ignore]
 elif model=="Recommendation":
  if collect_ignored:
   fields_to_ignore = {u"id", "creation_time"}
  else:
   fields_to_ignore = {u"id","user","lab_group", "saved", "model_version", "atoms", "creation_time", "nonsense"}
  dirty_fields = [field for field in Recommendation._meta.fields if field.name not in fields_to_ignore]
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

def get_ref_set(lab_group, reset_cache=True):
 ref_set = get_cache(lab_group, "DATAREFS")
 if not ref_set or reset_cache:
  ref_set = set(Data.objects.values_list('ref', flat=True))
  set_cache(lab_group, "DATAREFS", ref_set)
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
