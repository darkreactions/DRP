import subprocess
import uuid
import load_data
import rxn_calculator

POSITIVE = "4:4"
PREFIX = "/home/drp/web/darkreactions.haverford.edu/app/DRP/tmp/"
MODEL_BASE_DIR = "/home/drp/web/darkreactions.haverford.edu/app/DRP/models"
def gen_model(model_location):
	''' 
	gen_model("5.8.2014.UUID.model")
	'''

	rows = load_data.get_feature_vectors()
	headers = get_arff_headers()
	name = str(uuid.uuid4())
	print PREFIX + name + ".arff"
	with open(PREFIX+name + ".arff", "w") as raw:
		raw.write(headers+"\n")
		for row in rows:
			row[-1] = max(1, row[-1])
			raw.write(",".join([str(c) for c in row]) + "\n")
	
	subprocess.check_output("sh make_model.sh {0} {1}".format(MODEL_BASE_DIR + model_location, PREFIX+name+".arff"), shell=True)

def evaluate_model(results_location):
	with open(results_location) as results_file:
		for i in range(5):
			results_file.next()
		total = 0
		incorrect = 0
		false_positive = 0
		for row in results_file:
			if "\n" == row:
				continue
			row += 1
			if "+" in row:
				incorrect += 1
				if row.split()[2] == POSITIVE:
					false_positive += 1
	return incorrect, total, false_positive


def make_predictions(target_file, model_location):
	results_location = "/home/drp/web/darkreactions.haverford.edu/app/DRP/tmp/" + str(uuid.uuid4()) + ".out"
	subprocess.check_output("sh make_predictions.sh {0} {1}".format(target_file, model_location, results_location), shell=True)
	return results_location

def evaluate_real_results(lab_group = None):
	import sys, os
	sys.path.append('/home/drp/web/darkreactions.haverford.edu/app/DRP')
	os.environ.setdefault('DJANGO_SETTINGS_MODULE', 'DRP.settings')
	import DRP.models as m

	if lab_group:
		recommended = [o.outcome for o in m.Data.objects.filter(recommended="Yes", lab_group=lab_group)]
		unrecommended =  [o.outcome for o in m.Data.objects.filter(recommended="No", lab_group = lab_group)]
	else:
		recommended = [o.outcome for o in m.Data.objects.filter(recommended="Yes")]
		unrecommended =  [o.outcome for o in m.Data.objects.filter(recommended="No")]

	recommended_results = {'correct': recommended.count('4') + recommended.count('3'), 'incorrect': recommended.count('2') + recommended.count('1') + recommended.count('0') }
	unrecommended_results =  {'correct': unrecommended.count('4') + unrecommended.count('3'), 'incorrect': unrecommended.count('2') + unrecommended.count('1') + unrecommended.count('0') }
	return recommnded_results, unrecommened_results


def get_arff_headers():
	hdrs = rxn_calculator.headers
	XXX = 0
	for header in hdrs:
            if header[0:3] == "XXX":
                 XXX += 1
        headers = hdrs[XXX:]
	res = preface(headers, True)
	return res
	

def gen_specials():
	specials = {"outcome": "{1,2,3,4}", "slowCool": "{yes,no}",
		"leak": "{yes,no}", "purity": "{1,2}"}
	import rxn_calculator
	for bool_field in rxn_calculator.atomsz + rxn_calculator.bools:
		specials[bool_field] = "{yes,no}"
	return specials


def preface(headers, outcome = True, prefix = ""):
    res = "%  COMMENT \n%  NAME, DATE\n@relation rec_system" + prefix 
    specials = gen_specials()
    for header in headers:
        if header in specials.keys():
            if not (outcome and header == "purity") and not (not outcome and header == "outcome"):
                res += "\n@ATTRIBUTE " + header + " " + specials[header]
        else:
            res += "\n@ATTRIBUTE " + header + " NUMERIC"
    res += "\n\n@DATA\n"
    return res

