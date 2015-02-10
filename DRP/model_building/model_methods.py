import subprocess
import uuid
import rxn_calculator

import os, sys
full_path = os.path.dirname(os.path.realpath(__file__))+"/"
django_path = full_path[:full_path.rfind("/DRP/")]
if django_path not in sys.path:
  sys.path = [django_path] + sys.path
  os.environ['DJANGO_SETTINGS_MODULE'] = 'DRP.settings'

from DRP.retrievalFunctions import get_valid_data
from DRP.settings import MODEL_DIR, TMP_DIR

POSITIVE = "4:4"
INCORRECT = "+"

def makeBool(entry):
  if entry=="yes":
    return 1
  elif entry=="no":
    return 0
  else:
    return entry


def gen_model(model_name, description, data=None, clock=True, active=True, debug=False):
  '''
  gen_model("5.8.2014.model", "Some description of the model version.")
  will generate a model as the file "5.8.2014.model" and store the
  model's statistics in a ModelStats database entry.

  Optionally, only pre-specified data can be used to construct the model.
  '''

  from DRP.database_construction import store_ModelStats
  from DRP.fileFunctions import createDirIfNecessary
  from DRP.model_building.rxn_calculator import headers
  import time

  # Variable Setup
  correct_vals = None # `None` defaults to ["4:4","3:3"]
  name = str(int(time.time()))


  if not model_name or not description:
    raise Exception("Model needs a valid model_name and description!")


  #Make sure the model_name has no spaces in it.
  model_name = model_name.replace(" ","_").replace("\t","_").replace("/","-")
  print "Starting gen_model pipeline: {} ({})".format(model_name, name)


  if clock:
    import datetime
    print "Called gen_model at {}".format(datetime.datetime.now())


  # Get the valid reactions across all lab groups.
  print "Loading data entries."
  if data==None:
    data = get_valid_data()
    if debug:
      data = data[:100]
  print "... Loaded {} entries!".format(len(data))

  # Choose "training" and "test" data and construct a "sample model."
  #   From that sample model, see how well the actual model will perform.
  print "Creating a sample model to evaluate the model stats..."
  conf_json = sample_model_quality(data, name, clock=clock, debug=debug)

  print "Constructing the final ARFF file..."
  make_arff(name+"_final", data, clock=clock)

  #Using ALL of the data, now construct the full model.
  print "Building the actual model..."
  model_file_location = MODEL_DIR + model_name + ".model"
  arffFullName = TMP_DIR+name+"_final.arff"

  createDirIfNecessary(TMP_DIR)
  createDirIfNecessary(MODEL_DIR)

  tStart = time.clock()

  move = "cd {};".format(django_path)
  command = "bash DRP/model_building/make_model.sh"
  args = " {} {}".format(model_file_location, arffFullName)
  print "SUBPROCESS:\n{}".format(move+command+args)
  subprocess.check_output(move+command+args, shell=True)


  if clock: print "Took {} minutes.".format((time.clock()-tStart)/60.0)


  #Prepare a ModelStats entry and store it in the database.
  print "Creating a ModelStats entry in the database..."
  store_ModelStats(model_name, description, model_file_location, headers,
                   conf_json, correct_vals=correct_vals, active=active)

  print "Model generation successful..."




def get_current_model():
  from DRP.retrievalFunctions import get_latest_ModelStats
  model = get_latest_ModelStats()
  return model.filename


def map_to_zero_one(v):
	if v < 3:
		return 1
	else:
		return 2

# Builds and evaluates a sample model given some input data to be partitioned
#   into training and test data. Returns the stats of this sample model.
def sample_model_quality(data, name, clock=False, debug=False):
  from DRP.model_building.test_train_split import create_test_and_train_lists
  from DRP.model_building.load_data import create_reactant_keys
  from DRP.fileFunctions import createDirIfNecessary

  # Create reactant-combination keys for each data entry.
  dataKeys = create_reactant_keys(data)

  # Partitions the data/keys into separate test/training datasets.
  test, train = create_test_and_train_lists(data, dataKeys)

  # Generate the "train" and "test" ARFFs and give them a model-unique name.
  #TODO: ORIGINAL:   rows = [r[:-1] + [ map_to_zero_one(r[-1]) ] for r in rows]
  #TODO: Incorporate map_to_zero_one??
  #TODO: Why? Isn't used in FINAL make_arff?
  make_arff(name + "_test", test, clock=clock, debug=debug)
  make_arff(name + "_train", train, clock=clock, debug=debug)

  #Give the temporary ARFF files and the model respectable names.
  tmpPrefix = TMP_DIR + name
  modelFullName = MODEL_DIR + name + "_SAMPLE_MODEL"

  createDirIfNecessary(TMP_DIR)
  createDirIfNecessary(MODEL_DIR)

  # Start a new process to actually construct the model from the training data.
  move = "cd {};".format(django_path)
  command = "bash DRP/model_building/make_model.sh"
  args = " {} {}".format(modelFullName, tmpPrefix+"_train.arff")
  if debug:
    print "SUBPROCESS:\n{}".format(move+command+args)
  subprocess.check_output(move+command+args, shell=True)

  # Use the test data to make samplePredictions that can gauge the model quality.
  samplePredictions = make_predictions(tmpPrefix+"_test.arff", modelFullName, debug=debug)

  #Now that the sample model is created, evaluate its performance.
  return evaluate_results(samplePredictions)


def generate_weka_model(model_name, model_filepath):
  from DRP.settings import TMP_DIR
  import subprocess

  # Start a new process to actually construct the model from the training data.
  move = "cd {};".format(django_path)
  command = "bash DRP/model_building/make_model.sh"
  arff_filepath = TMP_DIR + model_name + "_train.arff"

  args = " {} {}".format(model_filepath, arff_filepath)
  subprocess.check_output(move+command+args, shell=True)



def evaluate_model(rows,keys):

	rows = [r[:-1] + [ map_to_zero_one(r[-1]) ] for r in rows]
	test,train = splitter.create_test_and_train_lists(rows, keys)


	name = str(uuid.uuid4())
	make_arff(name + "_test", test, True)
	make_arff(name + "_train", train, True)


	move = "cd {};".format(django_path)
	command = "bash DRP/model_building/make_model.sh"
	args = " {} {}".format(MODEL_DIR + name , TMP_DIR + name + "_train" + ".arff")
	subprocess.check_output(move+command+args, shell=True)
	results = make_predictions(TMP_DIR + name + "_test.arff", MODEL_DIR + name)

	performance, falsePositiveRate = evaluate_results(results)
	return performance, falsePositiveRate




def evaluate_results(results_location):
  """
  Sample structure of a results (tmp/*.out) file:

    inst#     actual     predicted     error    prediction
                         ...
    16        4:4        3:3           +        0.5
    17        3:3        3:3                    0.5
    18        3:3        3:3                    0.5
    19        3:3        3:3                    0.5
    20        1:1        3:3           +        0.5
    21        3:3        3:3                    0.5
    22        3:3        3:3                    0.5
    23        3:3        3:3                    0.5
    24        4:4        3:3           +        0.5
    25        4:4        4:4                    0.5
    26        3:3        3:3                    0.5
    27        2:2        4:4           +        0.5
    28        2:2        4:4           +        0.5
    29        4:4        4:4                    0.5

  """
  with open(results_location) as results_file:
    for i in range(5):
      results_file.next()

    conf_json = {}

    for row in results_file:
      if "\n" == row:
        continue

      actual = row.split()[1]
      guess = row.split()[2]

      # Create the confusion-dictionary...
      if actual in conf_json:
        if guess in conf_json[actual]:
          conf_json[actual][guess] += 1
        else:
          conf_json[actual][guess] = 1
      else:
        conf_json[actual] = {guess:1}

  return conf_json


def remove_indexes(row, unused):
  return [entry for i, entry in enumerate(row) if i not in unused]

def dict_to_list(calcDict, listFields):
  return [calcDict[field] for field in listFields]

# Creates ARFF contents given data and writes the contents to a file ('name').
def make_arff(name, data, clock=False, raw_list_input=False, debug=True, response="outcome"):
  import time
  from DRP.fileFunctions import createDirIfNecessary

  if clock and debug: tStart = time.clock()

  #Count the number of failed data (ie: invalid) that were passed to the make_arff.
  failed = 0
  i = 0

  createDirIfNecessary(TMP_DIR)
  fullFileName = TMP_DIR+name+".arff"

  if debug: print "Constructing: {}".format(fullFileName)

  with open(fullFileName, "w") as f:
    #Write the file headers.
    arff_fields, unused_indexes = get_used_fields()
    full_headers = preface(arff_fields, True)
    f.write(full_headers + "\n")

    #Write each datum to the file if possible.
    for datum in data:
      try:
        if raw_list_input:
          row = datum
        else:
          calcDict = datum.get_calculations_dict()
          if response=="outcome":
            calcDict[response] = max(1, int(calcDict[response]))
          row = dict_to_list(calcDict, arff_fields)

        #Write the row to the ARFF file.
        f.write(",".join([str(entry) for entry in row]) + "\n")

      except Exception as e:
        if debug: print "FAILED: {}".format(e)
        failed += 1
        #If the calculations_list failed/was invalid, erase it.
        datum.calculations = None
        datum.save()
      i += 1

  if debug: print "Completed: {} of {} data not usable.".format(failed, len(data))
  if clock and debug: print "Took {} seconds.".format(time.clock()-tStart)

  return fullFileName


def make_predictions(target_file, model_location, debug=False):
  results_location = "out".join(target_file.rsplit("arff", 1)) # Results file will be *.arff --> *.out
  move = "cd {};".format(django_path)
  comm = "bash DRP/model_building/make_predictions.sh"
  args = " {} {} {}".format(target_file, model_location, results_location)

  if debug: print "SUBPROCESS: {}".format(move+comm+args)
  subprocess.check_output(move+comm+args, shell=True)
  return results_location


#TODO: Rewrite this.
def evaluate_real_results(lab_group = None):
  import DRP.models as m

  if lab_group:
    recommended = [o.outcome for o in m.Data.objects.filter(recommended="Yes", lab_group=lab_group)]
    unrecommended =  [o.outcome for o in m.Data.objects.filter(recommended="No", lab_group = lab_group)]
  else:
    recommended = [o.outcome for o in m.Data.objects.filter(recommended="Yes")]
    unrecommended =  [o.outcome for o in m.Data.objects.filter(recommended="No")]

  recommended_results = {'correct': recommended.count('4') + recommended.count('3'), 'incorrect': recommended.count('2') + recommended.count('1') + recommended.count('0') }
  unrecommended_results =  {'correct': unrecommended.count('4') + unrecommended.count('3'), 'incorrect': unrecommended.count('2') + unrecommended.count('1') + unrecommended.count('0') }
  return recommended_results, unrecommended_results

def get_used_fields():
  unused_indexes = set()
  used = []
  for i, header in enumerate(rxn_calculator.headers):
    if header[0:3] == "XXX": # Don't use fields marked with an "XXX"
      unused_indexes.add(i)
    else:
      used.append(header)
  return used, unused_indexes


def removeUnused(row, unused_indexes=None):
  if unused_indexes is None:
    arff_fields, unused_indexes = get_used_fields()
  return [row[i] for i in xrange(len(row)) if i not in unused_indexes]


def get_good_result_tuples(results_location, rows, debug=False):
  reactions = []
  total = 0
  with open(results_location, "r") as results_file:
    # Remove the headers.
    for i in range(5):
      results_file.next()

    for row in results_file:
      if "+" not in row and row != "\n":
        clean = filter(lambda x:x!="" and x!="\n", row.split(" "))
        conf = float(clean[-1])
        index = int(clean[0])-1 # WEKA is 1-based, not 0-based.
        reactions.append((conf, rows[index]))
      total += 1

  if debug:
    "{} of {} reactions are good.".format(len(reactions), total)

  return reactions

#TODO:
"""
def get_arff_headers(zero_one=False):
	hdrs = rxn_calculator.headers
	XXX = 0
	for header in hdrs:
            if header[0:3] == "XXX":
                 XXX += 1
        headers = hdrs[XXX:]
        print "HEADERS: {}".format(len(headers))
	return res
"""

def gen_specials():
	import rxn_calculator

	specials = {"outcome": "{1,2,3,4}", "slowCool": "{yes,no}",
		"leak": "{yes,no}"}

	for bool_field in rxn_calculator.atomsz + rxn_calculator.bools:
		specials[bool_field] = "{yes,no}"

	return specials


def preface(headers, data):
    res = "%  COMMENT \n%  NAME, DATE\n@relation rec_system"
    specials = gen_specials()
    for header in headers:
        if header in specials.keys():
            res += "\n@ATTRIBUTE " + header + " " + specials[header]
        else:
            res += "\n@ATTRIBUTE " + header + " NUMERIC"
    res += "\n\n@DATA\n"
    return res

