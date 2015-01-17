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
  from DRP.chemspider import search_chemspider

  try:
    if not entry.custom: #Only update compounds that are not custom.
      #Get the most up-to-date ChemSpider info for a given CAS/compound.
      query = search_chemspider([entry.CAS_ID, entry.compound])
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
  reaction.get_atoms(refresh=True)

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
    clean_data["image_url"], clean_data["smiles"], clean_data["mw"] = query.imageurl, query.smiles, query.molecularweight

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

