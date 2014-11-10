from retrievalFunctions import *
from logPrinting import print_error, print_log

   # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
   # # # # # # # # # # # # #  CG_calculations  # # # # # # # # # # # # # # # #
   # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


   # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
   # # # # # # # # # # # # #  Recommendations  # # # # # # # # # # # # # # # #
   # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

def store_ModelStats(model_name, model_description, filename, true_pos, true_neg, false_pos, false_neg, active=True):

  from DRP.models import ModelStats
  import datetime

  #Create an instance of the model.
  model = ModelStats()
  model.datetime = datetime.datetime.now()

  model.set_values(false_neg, false_pos, true_neg, true_pos)
  model.title = model_name
  model.description = model_description
  model.filename = filename
  model.active = active

  #Save the model.
  model.save()

  return model

#Organize the reactants into sublists: [[reactant, quantity, unit], ...]
def partition_reactant_fields(lst):
 num_fields = CONFIG.fields_per_reactant
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
def get_vars_from_list(lst):
 verbose_headers = get_model_field_names(verbose=True, model="Data")[16:]
 info = [[i, j] for (i,j) in zip(verbose_headers, lst)]
 return info

#Creates the Recommendation entry, but does not store it in database.
def field_list_to_Recommendation(lab_group, lst, in_bulk=False, debug=False):
 from DRP.models import Recommendation, get_model_field_names, get_atom_set_from_abbrevs, get_abbrevs_from_reaction
 import datetime

 try:
  new_rec = Recommendation()
  #Set the self-assigning fields:
  setattr(new_rec, "lab_group", lab_group)
  setattr(new_rec, "score", float(lst[0]))

  #Set the non-user field values.
  fields = get_model_field_names(model="Recommendation")

  for (field, value) in zip(fields, lst[2:]): #Ignore the reference field.

   #Translate Booleans into Boolean values.
   if field in bool_fields:
    value = True if value[0].lower() in "1tyc" else False

   if debug:
     print "... Setting '{}' as '{}' ({})".format(field, value, type(value))

   setattr(new_rec, field, value)

  new_rec.atoms = "".join(get_atom_set_from_abbrevs(lab_group, get_abbrevs_from_reaction(new_rec)))

  if not in_bulk:
    new_rec.date_dt = datetime.datetime.now()

  return new_rec
 except Exception as e:
  raise Exception("Recommendation construction failed: {}".format(e))

def store_new_RankedReaction(json_input):
 #Allow either valid json/dicts or strings that can be loaded as such.
 if type(json_input)==str:
  json_input = json.loads(json_input)

 try:
  rxn_list = RankedReactionList()
  rxn_list.seed = json.dumps(json_input["seed"])
  rxn_list.original_list = json.dumps(json_input["targets"])
  rxn_list.save()
  return rxn_list
 except Exception as e:
  raise Exception("RankedReactionList construction failed: {}".format(e))

def store_new_RankedReaction_list(list_of_rankedrxn_lists):
 #Variable Setup:
 successes = 0
 count = 0

 for reaction_json in list_of_rankedrxn_lists:
  try:
   count += 1
   store_new_RankedReaction(reaction_json)
   successes += 1
  except Exception as e:
    print "RankedReactionList {} could not be constructed: {}".format(count, e)

 print "Finished creating and storing {} of {} items!.".format(successes, count)


def store_new_Recommendation_list(lab_group, recommendations, version_notes = "", seed_source=None, debug=True):
  from DRP.models import get_Lab_Group, ModelStats
  from DRP.retrievalFunctions import get_latest_ModelStats
  import datetime

  lab_group = get_Lab_Group(lab_group)

  call_time = datetime.datetime.now()

  #Either prepare a new "Model Version" or use the latest one for a lab group.
  try:
    model = get_latest_ModelStats(lab_group)
  except Exception as e:
    print_error("Model not gathered for the Recommendation list: {}".format(e))

  #Store the actual Recommendation entries.
  num_success = 0
  count = 0
  for rec in recommendations:
    count += 1
    try:
     print "1: {}".format(rec)
     new_rec = field_list_to_Recommendation(lab_group, rec, in_bulk=True)
     new_rec.date_dt = call_time
     new_rec.model_version = model
     if seed_source:
       new_rec.seeded = True
       new_rec.seed = seed_source #Record if this recommendation is seeded.

     new_rec.save() #Store this recommendation in the database
     num_success += 1
    except Exception as e:
      print_error("Recommendation {} could not be constructed: {}".format(count, e))

  if debug: print_log("Finished creating and storing {} of {} items!.".format(num_success, count))
