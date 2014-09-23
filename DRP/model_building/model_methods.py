import subprocess
import uuid
import time
import load_data
import rxn_calculator

import sys, os
django_dir = os.path.dirname(os.path.realpath(__file__)).split("DRP")[0]
django_path = "{}/DRP".format(django_dir)
if django_path not in sys.path:
  sys.path.append("{}/DRP".format(django_dir))

os.environ['DJANGO_SETTINGS_MODULE'] = 'DRP.settings'

from DRP.database_construction import store_ModelStats
from DRP.retrievalFunctions import get_valid_data
from DRP.settings import BASE_DIR, MODEL_DIR, TMP_DIR

POSITIVE = "2:2"

def gen_model(model_name, description, data=None):
  '''
  gen_model("5.8.2014.model", "Some description of the model version.")
  will generate a model as the file "5.8.2014.model" and store the
  model statistics in a ModelStats database entry.

  Optionally, only certain data will be used to construct the model.
  '''

  if not model_name or not description:
    raise Exception("Model needs a valid model_name and description!")

  #Make sure the model_name has no spaces in it.
  model_name = model_name.replace(" ","_")
  name = str(int(time.time()))
  print "Constructing model: {} ({})".format(model_name, name)

  # Get the valid reactions across all lab groups.
  print "Loading data entries."
  if not data:
    data = get_valid_data()

  # Choose "training" and "test" data and construct a "sample model."
  #   From that sample model, see how well the actual model will perform.
  print "Creating a sample model to evaluate the model stats..."
  truePositiveRate, falsePositiveRate = sample_model_quality(data, name)

  print "Constructing the final ARFF file..."
  make_arff(name+"_final", data)

  #Using ALL of the data, now construct the full model.
  print "Building the actual model..."
  modelFullName = MODEL_DIR + model_name
  arffFullName = TMP_DIR+name+"_final.arff"
  command = "bash DRP/model_building/make_model.sh"
  args = " {} {}".format(modelFullName, arffFullName)
  print "COMMAND: {} {}".format(command, args)
  subprocess.check_output(command+args, shell=True)


  #Prepare a ModelStats entry and store it in the database.
  print "Creating a ModelStats entry in the database..."
  update_dashboard(false_positive = falsePositiveRate,
                   model_performance = truePositiveRate,
                   description=description,
                   model_name = model_name)

  print "Model generation successful..."




def get_current_model():
	return MODEL_DIR + sorted([f for f in os.listdir(MODEL_DIR) if "model" in f], key = lambda x: x.split(".")[0], reverse = True)[0]


def map_to_zero_one(v):
	if v < 3:
		return 1
	else:
		return 2

# Builds and evaluates a sample model given some input data to be partitioned
#   into training and test data. Returns the stats of this sample model.
from DRP.model_building.test_train_split import create_test_and_train_lists
from DRP.model_building.load_data import create_reactant_keys
def sample_model_quality(data, name):

  # Create reactant-combination keys for each data entry.
  dataKeys = create_reactant_keys(data)

  # Partitions the data/keys into separate test/training datasets.
  test, train = create_test_and_train_lists(data, dataKeys)

  # Generate the "train" and "test" ARFFs and give them a model-unique name.
  #TODO: ORIGINAL:   rows = [r[:-1] + [ map_to_zero_one(r[-1]) ] for r in rows]
  #TODO: Incorporate map_to_zero_one??
  #TODO: Why? Isn't used in FINAL make_arff?
  make_arff(name + "_test", test)
  make_arff(name + "_train", train)

  #Give the temporary ARFF files and the model respectable names.
  tmpPrefix = TMP_DIR + name
  fullModelName = MODEL_DIR + name + "_SAMPLE_MODEL"

  # Start a new process to actually construct the model from the training data.
  command = "bash DRP/model_building/make_model.sh"
  args = " {} {}".format(fullModelName,  tmpPrefix+"_train.arff")
  print "COMMAND: {} {}".format(command, args)
  subprocess.check_output(command+args, shell=True)

  # Use the test data to make samplePredictions that can gauge the model quality.
  samplePredictions = make_predictions(tmpPrefix+"_test.arff", fullModelName)

  #Now that the sample model is created, evaluate its performance.
  truePositiveRate, falsePositiveRate = evaluate_results(samplePredictions)
  return truePositiveRate, falsePositiveRate




def evaluate_model(rows,keys):

	rows = [r[:-1] + [ map_to_zero_one(r[-1]) ] for r in rows]
	test,train = splitter.create_test_and_train_lists(rows, keys)


	name = str(uuid.uuid4())
	make_arff(name + "_test", test, True)
	make_arff(name + "_train", train, True)


	subprocess.check_output("bash DRP/model_building/make_model.sh {0} {1}".format(MODEL_DIR + name , TMP_DIR + name + "_train" + ".arff"), shell=True)
	results = make_predictions(TMP_DIR + name + "_test.arff", MODEL_DIR + name)

	performance, falsePositiveRate = evaluate_results(results)
	return performance, falsePositiveRate




def evaluate_results(results_location):
	with open(results_location) as results_file:
		for i in range(5):
			results_file.next()
		total = 0
		incorrect = 0
		false_positive = 0
		true_negative = 0
		for row in results_file:
			if "\n" == row:
				continue
			total += 1
			if "+" in row: #If it's "unrecommended"...
				incorrect += 1
				if row.split()[2] == POSITIVE:
					false_positive += 1
				else:
					true_negative += 1

	#Calculate the rates from the various counts.
	falsePositiveRate = false_positive/float(true_negative+false_positive) if (true_negative + false_positive) else 0
	truePositiveRate = (total - incorrect) / float(total) if total else 0

	return truePositiveRate, falsePositiveRate



# Creates ARFF contents given data and writes the contents to a file ('name').
def make_arff(name, data):
  #Get the valid headers.
  headers = get_arff_headers()

  print "ARFFING"

  #Count the number of failed data (ie: invalid) that were passed to the make_arff.
  failed = 0
  i = 0

  fullFileName = TMP_DIR+name+".arff"
  with open(fullFileName, "w") as f:
    #Write the file headers.
    f.write(headers + "\n")
    print len(data),
    #Write each datum to the file if possible.
    for datum in data:
      print i
      try:
        row = datum.get_calculations_list()

        #TODO: Figure out a better way to get this into the form Paul's
        #      scripts need...
        # Remove the "raw" data before ARFFing (that is, the data
        #   from the "ref" up to the first calculated field).
        row = row[19:-2] + [row[-1]]
        #Write the row to the ARFF file.
        row[-1] = max(1, row[-1])
        f.write(",".join([str(entry) for entry in row]) + "\n")

      except Exception as e:
        print e
        failed += 1
        #If the calculations_list failed/was invalid, erase it.
        datum.calculations = None
        datum.save()
      i += 1

  print "{} of {} data not usable by make_arff".format(failed, len(data))


#TODO: Delete
"""
def make_arff(name, rows, zero_one = False):
	headers = get_arff_headers(zero_one)

	with open(TMP_DIR+name + ".arff", "w") as raw:
		raw.write(headers+"\n")
		for row in rows:
			row[-1] = max(1, row[-1])
			raw.write(",".join([str(c) for c in row]) + "\n")
"""



def make_predictions(target_file, model_location):
  results_location = TMP_DIR + str(uuid.uuid4()) + ".out"
  comm = "bash DRP/model_building/make_predictions.sh {0} {1} {2}".format(target_file, model_location, results_location)
  subprocess.check_output(comm, shell=True)
  return results_location


def update_dashboard(
  false_positive = None,
  model_performance = None,
  rec_estimate = None,
  empirical_success = None,
  description = "",
  model_name = ""):

  import DRP.models as m

  #If a specific datum is missing for some reason, use the previous one.
  last = m.ModelStats.objects.last()
  if false_positive is None:
    false_positve = last.false_positive_rate
  if model_performance is None:
    model_performance = last.performance
  if rec_estimate is None:
    rec_estimate = last.estimated_success_rate
  if empirical_success is None:
    empirical_success, _ = evaluate_real_results()
    empirical_success = empirical_success['correct'] / float(empirical_success['correct'] + empirical_success['incorrect'])

  if not model_name:
    model_name = "Automated update"

        #Store these stats in the database
  store_ModelStats(false_positive, empirical_success, rec_estimate, model_performance, description, model_name)


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


def get_arff_headers(zero_one=False):
	hdrs = rxn_calculator.headers
	XXX = 0
	for header in hdrs:
            if header[0:3] == "XXX":
                 XXX += 1
        headers = hdrs[XXX:]
	res = preface(headers, True, zero_one = zero_one)
	return res


def gen_specials(zero_one=False):
	specials = {"outcome": "{1,2,3,4}", "slowCool": "{yes,no}",
		"leak": "{yes,no}", "purity": "{1,2}"}

	if zero_one:
		specials["outcome"] = "{1,2}"
	import rxn_calculator
	for bool_field in rxn_calculator.atomsz + rxn_calculator.bools:
		specials[bool_field] = "{yes,no}"
	return specials


def preface(headers, outcome = True, prefix = "", zero_one = False):
    res = "%  COMMENT \n%  NAME, DATE\n@relation rec_system" + prefix
    specials = gen_specials(zero_one)
    for header in headers:
        if header in specials.keys():
            if not (outcome and header == "purity") and not (not outcome and header == "outcome"):
                res += "\n@ATTRIBUTE " + header + " " + specials[header]
        else:
            res += "\n@ATTRIBUTE " + header + " NUMERIC"
    res += "\n\n@DATA\n"
    return res

