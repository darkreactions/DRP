from retrievalFunctions import *
from DRP.logPrinting import print_error, print_log

   # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
   # # # # # # # # # # # # #  CG_calculations  # # # # # # # # # # # # # # # #
   # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


   # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
   # # # # # # # # # # # # #  Recommendations  # # # # # # # # # # # # # # # #
   # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

def store_ModelStats(falsePositive, actualSuccess, estimatedSuccess, modelPerformance, description, datetime_input=None):
  #Create an instance of the model.
  model_stats = ModelStats()

  #Set the fields of model_stats object.  
  model_stats.false_positive_rate = falsePositive
  model_stats.actual_success_rate = actualSuccess
  model_stats.estimated_success_rate = estimatedSuccess
  model_stats.performance = modelPerformance
  model_stats.description = description

  if (not datetime_input):
    datetime_input = datetime.datetime.now()


  model_stats.datetime = datetime_input

  #Save the model.
  model_stats.save()  

  return model_stats

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
def field_list_to_Recommendation(lab_group, lst, in_bulk=False):
 try:
  new_rec = Recommendation()
  #Set the self-assigning fields:
  setattr(new_rec, "lab_group", lab_group)
  setattr(new_rec, "score", float(lst[0]))

  #Set the non-user field values.
  fields = get_model_field_names(model="Recommendation")
  for (field, value) in zip(fields, lst[2][1:]): #Ignore the reference field.
   #Translate Booleans into Boolean values.
   if field in bool_fields:
    value = True if value[0].lower() in "1tyc" else False
   setattr(new_rec, field, value)
   
  new_rec.atoms = "".join(get_atom_set_from_abbrevs(lab_group, get_abbrevs_from_reaction(new_rec)))
 
  if not in_bulk:
    new_rec.date = str(datetime.datetime.now())

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


def store_new_Recommendation_list(lab_group, list_of_recommendations, version_notes = "", seed_source=None, new_model=False):
  lab_group = get_Lab_Group(lab_group)
  
  call_time = str(datetime.datetime.now())
 
  #Either prepare a new "Model Version" or use the latest one for a lab group. 
  try:
    if new_model:
      model = Model_Version()
      model.model_type = "Recommendation"
      model.date = call_time
      model.notes = version_notes
      model.lab_group = lab_group
      model.save()
    else:
      model = get_latest_Model_Version(lab_group)
  except Exception as e:
    print_error("Model not gathered for the Recommendation list: {}".format(e))

  #Store the actual Recommendation entries.
  num_success = 0
  count = 0
  for i in list_of_recommendations:
    count += 1
    try:
     new_rec = field_list_to_Recommendation(lab_group, i, in_bulk=True)
     new_rec.date = call_time
     new_rec.model_version = model
     if seed_source:
       new_rec.seeded = True
       new_rec.seed = seed_source #Record if this recommendation is seeded.
  
     new_rec.save() #Store this recommendation in the database
     num_success += 1
    except Exception as e:
      print_error("Recommendation {} could not be constructed: {}".format(count, e))
  
  print_log("Finished creating and storing {} of {} items!.".format(num_success, count))
