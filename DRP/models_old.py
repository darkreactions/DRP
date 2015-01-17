from django.forms import *
from django.core import validators
from django.contrib.auth.models import User, Group
from django.db import models
from django.db.models import Q

from data_config import CONFIG
from validation import *
from uuid import uuid4
from CGCalculator import CGCalculator
from collections import defaultdict
from subprocess import Popen
from settings import LOG_DIR, BASE_DIR, MODEL_DIR
from cacheFunctions import get_cache, set_cache

import json, random, string, datetime, operator
import chemspipy



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
 import rdkit.Chem as Chem
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
 from DRP.data_config import CONFIG
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


def refresh_compound_guide(lab_group=None, verbose=False, debug=False, clear=False):
  #Either refresh all of the data or the data for a specific lab group.
  if lab_group:
    query = get_lab_CG(lab_group)
  else:
    query = CompoundEntry.objects.all()

  if clear:
    if debug: print "Clearing all CG_calculations..."
    CompoundEntry.objects.all().update(calculations=None)
    CompoundEntry.objects.all().update(calculations_failed=False)
    CG_calculations.objects.all().delete()

  #Actually perform the refresh/update.
  for i, compound in enumerate(query):
    try:
      if verbose and i%10==0: print "... {}.".format(i)
      update_compound(compound, debug=debug)
    except Exception as e:
      if debug: print "Could not update: {}\n\t".format(compound, e)


#Update the compound by reloading the ChemSpider search data.
def update_compound(entry, debug=False):
  from DRP.fileFunctions import createDirIfNecessary

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
        if debug:
          print "Found legacy entry (should be `custom`): {}".format(entry.compound)

    else:
        perform_calcs = False
        entry.calculations = None
        entry.image_url, entry.smiles, entry.mw = "","",""
    entry.save()

    if perform_calcs:

      #Start a new compound-calc worker to determine compound properties.
      comp_log_dir = LOG_DIR+"/compound_calculations"
      createDirIfNecessary(comp_log_dir)

      err_log = open(comp_log_dir+"/error.log","a")
      act_log = open(comp_log_dir+"/process.log","a")
      worker_script = BASE_DIR+"/DRP/compound_calculations/calculate_compound_properties.py"
      command = "python {} {}".format(worker_script, entry.id)
      #Log to the files above and make the worker independent of the parent process.
      Popen(command.split(), stdout=act_log, stderr=err_log, close_fds=True)

  except Exception as e:
    entry.calculations_failed = False
    entry.save()
    raise Exception("Compound update ({}) failed: {}".format(entry, e))

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


############### USER and LAB INTEGRATION #######################
############### USER CREATION #######################
############### DATA ENTRY ########################


# Convert any number-like strings to floats.
def make_float(string):
  try:
    return float(string)
  except:
    return string

############### RECOMMENDATIONS ########################
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

def create_CG_calcs_if_needed(compound, smiles, compound_type, ):
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

def perform_CG_calculations(only_missing=True, lab_group=None, attempt_failed = True, verbose=False):
 #Variable Setup
 success = 0
 i = 0

 cg = CompoundEntry.objects.all()
 if only_missing:
  cg = cg.filter(calculations=None)
 if lab_group:
  cg = cg.filter(lab_group=lab_group)
 if not attempt_failed:
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

def get_Data_with_compound(compound):
 Q_list = [Q(("reactant_{}".format(i),compound)) for i in CONFIG.reactant_range()]
 return Data.objects.filter(reduce(operator.or_, Q_list))


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
