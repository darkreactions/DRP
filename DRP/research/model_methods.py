import subprocess
import uuid

POSITIVE = "4:4"
def gen_model(train_arff, model_location):
	subprocess.check_output("sh make_model.sh {0} {1}".format(model_location, train_arff), shell=True)	

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
