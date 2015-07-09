from retrievalFunctions import *
from DRP.validation import bool_fields
from logPrinting import print_error, print_log
from DRP.retrievalFunctions import get_compound_by_name
import sys

   # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
   # # # # # # # # # # # # #  CG_calculations  # # # # # # # # # # # # # # # #
   # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


   # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
   # # # # # # # # # # # # #  Recommendations  # # # # # # # # # # # # # # # #
   # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

def store_ModelStats(model_name, model_description, filename, model_fields,
                     conf_json, correct_vals=None, active=True):

  from models import ModelStats
  import datetime

  #Create an instance of the model.
  model = ModelStats()
  model.datetime = datetime.datetime.now()

  model.set_confusion_table(conf_json)
  model.set_correct_vals(correct_vals)
  model.set_used_fields(model_fields)

  model.title = model_name
  model.description = model_description
  model.filename = filename
  model.active = active

  #Save the model.
  model.save()

  model.summary()

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
    from DRP.models import Recommendation, get_model_field_names
    import datetime

  # Try catch commented for debugging purposes
  #try:
    debug=True
    new_rec = Recommendation()
    #Set the self-assigning fields:
    setattr(new_rec, "lab_group", lab_group)
    setattr(new_rec, "score", float(lst[0]))

    #Set the non-user field values.
    fields = get_model_field_names(model="Recommendation")
    if debug: sys.stdout.write("fields: " + str(fields) + "\n:")
    if debug: sys.stdout.write("len(fields): " + str(len(fields)) + "\n:")
    if debug: sys.stdout.flush()

    for (field, value) in zip(fields[1:], lst[2:]): #Ignore the reference field (and conf)

      #Translate Booleans into Boolean values.
      # Assuming this was supposed to be the bool_fields in DRP.validation and adding import to top of file (there was no value assigned to the name bool_fields)
      if field in bool_fields:
        value = True if value[0].lower() in "1tyc" else False

      # Use the CompoundEntry object to build the rec instead of the compound string
      if field.startswith("reactant_fk"):
        value = get_compound_by_name(value)

      if debug: sys.stdout.write("jkl; graaaaaaaaaaaaaaaa \n:")
      if debug: sys.stdout.flush()

      if debug:
        print "... Setting '{}' as '{}' ({})".format(field, value, type(value))

      if debug: sys.stdout.write("jkl; gruuuuuuuuuuuuuuuu \n:")
      if debug: sys.stdout.flush()

      if debug: sys.stdout.write("new_rec, field, value: " + str(new_rec) + ", " + str(field) + ", " + str(value) + "\n")
      if debug: sys.stdout.flush()
      setattr(new_rec, field, value)

    if debug: sys.stdout.write("asdf new_rec.atoms: " + str("".join(new_rec.get_atoms())) + "\n:")
    if debug: sys.stdout.write("asdf type(new_rec.atoms): " + str(type("".join(new_rec.get_atoms()))) + "\n:")
    if debug: sys.stdout.flush()
    new_rec.atoms = "".join(new_rec.get_atoms())

    if not in_bulk:
      new_rec.date_dt = datetime.datetime.now()

    return new_rec
  #except Exception as e:
    import os
    os.system('echo "' + "debug point" + '"|espeak')

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

  # I believe the following comment is no longer applicable to post-ModelStats refactoring code...
  #Either prepare a new "Model Version" or use the latest one for a lab group.
  try:
    # Was `model = get_latest_ModelStats(lab_group)`, but get_latest_ModelStats takes no argument.
    # Removing argument. Should get_latest_ModelStats have an optional lab_group argument or was this a
    # refactoring mistake?
    model = get_latest_ModelStats()
  except Exception as e:
    print_error("Model not gathered for the Recommendation list: {}".format(e))

  #Store the actual Recommendation entries.
  num_success = 0
  count = 0
  for rec in recommendations:
      count += 1
    # Try-catch commented for debugging purposes
    #try:
      print "1: {}".format(rec)
      if debug: print "length: ", len(rec)
      new_rec = field_list_to_Recommendation(lab_group, rec, in_bulk=True)
      new_rec.date_dt = call_time
      new_rec.model_version = model
      if seed_source:
        new_rec.seeded = True
        new_rec.seed = seed_source #Record if this recommendation is seeded.

      new_rec.save() #Store this recommendation in the database
      num_success += 1
    #except Exception as e:
      print_error("Recommendation {} could not be constructed: {}".format(count, e))

  if debug: print_log("Finished creating and storing {} of {} items!.".format(num_success, count))
