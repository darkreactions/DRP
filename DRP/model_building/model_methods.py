import subprocess
import uuid
import load_data
import rxn_calculator

import sys, os
sys.path.append('/home/drp/web/darkreactions.haverford.edu/app/DRP')
os.environ.setdefault('DJANGO_SETTINGS_MODULE', 'DRP.settings')


import DRP.database_construction as dbc


POSITIVE = "2:2"
PREFIX = "/home/drp/web/darkreactions.haverford.edu/app/DRP/tmp/"
MODEL_BASE_DIR = "/home/drp/web/darkreactions.haverford.edu/app/DRP/models/"
def gen_model(model_name):
	''' 
	gen_model("5.8.2014.UUID.model")
	will generate a model in 5.8.2014.UUID.model
	'''

	name = str(uuid.uuid4())

	rows, keys = load_data.get_feature_vectors(keys=True)
	print "evaluating model"
	inc, total, false_p =  evaluate_model(rows, keys)
	make_arff(name, rows)
	
	subprocess.check_output("sh make_model.sh {0} {1}".format(MODEL_BASE_DIR + model_name, PREFIX+name+".arff"), shell=True)

	performance = (total - inc) / float(total) if total != 0 else 0

	update_dashboard(false_positive = false_p, model_performance = performance, model_name = model_name)

	


def get_current_model():
	return sorted([f for f in os.listdir("/home/drp/web/darkreactions.haverford.edu/app/DRP/models/") if "." in f], key = lambda x: x.split(".")[0], reverse = True)[0]
	

def map_to_zero_one(v):
	if v < 3:
		return 1 
	else:
		return 2


def evaluate_model(rows,keys):
	import test_train_split as splitter

	rows = [r[:-1] + [ map_to_zero_one(r[-1]) ] for r in rows]
	test,train = splitter.create_test_and_train_lists(rows, keys)


	name = str(uuid.uuid4())
	make_arff(name + "test", test, True)
	make_arff(name + "train", train, True)
	
	subprocess.check_output("sh make_model.sh {0} {1}".format(MODEL_BASE_DIR + name , PREFIX + name + "train" + ".arff"), shell=True)
	results = make_predictions(PREFIX + name + "test.arff", MODEL_BASE_DIR + name)
	incorrect, total, false_positive = evaluate_results(results) 
	
	#TODO: add new entry
	return incorrect, total, false_positive

	
	 

def evaluate_results(results_location):
	with open(results_location) as results_file:
		for i in range(5):
			results_file.next()
		total = 0
		incorrect = 0
		false_positive = 0
		for row in results_file:
			if "\n" == row:
				continue
			total += 1
			if "+" in row:
				incorrect += 1
				if row.split()[2] == POSITIVE:
					false_positive += 1
	return incorrect, total, false_positive


def make_arff(name, rows, zero_one = False):


	headers = get_arff_headers(zero_one)

	print PREFIX + name + ".arff"
	with open(PREFIX+name + ".arff", "w") as raw:
		raw.write(headers+"\n")
		for row in rows:
			row[-1] = max(1, row[-1])
			raw.write(",".join([str(c) for c in row]) + "\n")




def make_predictions(target_file, model_location):
	results_location = "/home/drp/web/darkreactions.haverford.edu/app/DRP/tmp/" + str(uuid.uuid4()) + ".out"
	subprocess.check_output("sh make_predictions.sh {0} {1} {2}".format(target_file, model_location, results_location), shell=True)
	return results_location


def update_dashboard(
	false_positive = None, 
	model_performance = None, 
	rec_estimate = None, 
	empirical_success = None,
	model_name = None): 

	import DRP.models as m

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

	if model_name is None:
		model_name = "Automated update"

	print false_positive, model_performance, rec_estimate, empirical_success, model_name
	dbc.store_ModelStats(false_positive, empirical_success, rec_estimate, model_performance, model_name) 
	
	
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

