import  json, uuid, subprocess
from DRP.research import metrics, load_cg, parse_rxn, clean2arff, rebuildCDT
import sys, os

sys.path.append('/home/drp/web/darkreactions.haverford.edu/app/DRP')
os.environ.setdefault('DJANGO_SETTINGS_MODULE', 'DRP.settings')
import DRP.models
from DRP.settings import TMP_DIR, BASE_DIR
from DRP.errorReporting import print_error

sim = metrics.Euclidean("euclidean", DRP.models).sim

#Variable Setup.
CG = load_cg.get_cg()

reactant_fields = {
	"reactant_1": ("quantity_1", "unit_1"),
	"reactant_2": ("quantity_2", "unit_2"),
	"reactant_3": ("quantity_3", "unit_3"),
	"reactant_4": ("quantity_4", "unit_4"),
	"reactant_5": ("quantity_5", "unit_5"),
}


pH_range = [1,3,5]
number_of_amine_mole_steps = 10
steps_per_amine = len(pH_range)*number_of_amine_mole_steps 

def get_reactants_indices(reaction):
	indices = dict()
	for r_f in reactant_fields:
		if getattr(reaction, r_f) not in CG:
			#print "skipping: {0}".format(getattr(reaction, r_f))
			continue

		reactant_type = is_valid_reactant(reaction, r_f)
		if reactant_type == "Org":
			indices["org"] = r_f 
		elif reactant_type == "Inorg":
			if "metal_1" in indices:
				indices["metal_2"] = r_f
			else:
				indices["metal_1"] = r_f 
		elif reactant_type == "Water":
			indices["water"] = r_f	
	if "metal_1" not in indices or "metal_2" not in indices or "org" not in indices or "water" not in indices:
		raise Exception("Missing a metal or amine: " + str(indices))
	return indices


def is_valid_reactant(reaction, r_f):
	if getattr(reaction, r_f) not in CG:
		raise Exception("Reactant not in cg: {0}".format(getattr(reaction,r_f)))
	if getattr(reaction, reactant_fields[r_f][1]) != "g":
		raise Exception("These aren't grams: {0}".format(getattr(reaction,r_f)))
	if not CG[getattr(reaction,r_f)]["mw"]:
		raise Exception("MW is zero: {0}".format(getattr(reaction,r_f)))

	return CG[getattr(reaction,r_f)]["type"]

def get_amine_moles(reaction, amine_index):
	return CG[getattr(reaction,amine_index)]["mw"] / float(getattr(reaction,reactant_fields[amine_index][0]))
	

def get_amine_range(moles):
	return [i/5.0 for i in range(11)]

def row_generator(reaction, indices, amine_moles, amine_list):
	
	metal_1 = getattr(reaction,indices["metal_1"])
	metal_2 = getattr(reaction,indices["metal_2"])
	amine = getattr(reaction,indices["org"])
	water = getattr(reaction,indices["water"])

	metal_1_mass = getattr(reaction,reactant_fields[indices["metal_1"]][0])
	metal_2_mass = getattr(reaction,reactant_fields[indices["metal_2"]][0])
	amine_moles = get_amine_moles(reaction, indices["org"])
	water_mass = getattr(reaction,reactant_fields[indices["water"]][0])

	amine_range = get_amine_range(amine_moles)
	pH_range = [1,3,5]

	for amine in amine_list:
		if amine not in CG:
			print "Not in CG: {0}".format(amine)
			continue
		for moles in amine_range:
			mass = moles*CG[amine]["mw"]
			for pH in pH_range:
				yield ["--", metal_1, metal_1_mass, "g", metal_2, metal_2_mass,
					"g", amine, mass, "g", water, water_mass, "g", "","","", 
					getattr(reaction,"temp"), getattr(reaction,"time"), pH, "yes", "no", 4, 2, ""]

def get_candidates(results, idx, raw_rows):
	candidates = []
	for i in range(steps_per_amine):
		if i + idx in results:
			candidates.append( (results[i+idx], raw_rows[i+idx]))
	if len(candidates) == 0:
		return None, 0.0
	best_conf = max(candidates, key=lambda x: x[0])[0]
	candidates = filter(lambda x: x[0] == best_conf, candidates)
	candidates = [c[1] for c in candidates]

	avg_mass = 0.0
	avg_ph = 0.0
	for r in candidates:
		avg_mass += float(r[8])
		avg_ph += float(r[18])
	
	avg_mass = avg_mass / float(len(candidates))
	avg_ph = avg_ph / float(len(candidates))

	smallest_dist = 10000000000000000000
	best_row = None

	for row in candidates:
		dist = (avg_mass-float(row[8]))**2 + (avg_ph - float(r[18]))**2
		if dist < smallest_dist:
			smallest_dist = dist
			best_row = row

	assert(best_row is not None)

	return best_row, best_conf

def generate_grid(reaction, amine_list, debug=True):
	#Variable Setup.
	indices = get_reactants_indices(reaction)
	amine_moles = get_amine_moles(reaction, indices["org"])
	prefix = TMP_DIR
	fileprefix = str(uuid.uuid4())
	ml_convert = json.load(open(BASE_DIR+"/DRP/research/mlConvert.json"))
	hdrs = ",".join(rebuildCDT.headers)

	row_gen = row_generator(reaction, indices, amine_moles, amine_list)
	raw_rows = [row for row in row_gen]
	#print raw_rows
	rows =[hdrs] + [ ",".join([str(c).replace(",","c") for c in parse_rxn.parse_rxn(row, CG, ml_convert)]) for row in raw_rows] 

	with open(prefix + fileprefix + ".csv", "w") as outfile:
		for row in rows:
			outfile.write(row+"\n")

	if debug: print prefix + fileprefix + ".csv"
			
	clean2arff.clean(prefix+fileprefix)

	#TODO: rewrite test_model
	cmd = "sh {0}/DRP/research/test_model.sh {1}".format(BASE_DIR, fileprefix)
	result = subprocess.check_output(cmd, shell=True)

	if debug: print result, cmd

	weka_results_file = fileprefix + ".out" 
	results = dict()
	with open(prefix + weka_results_file, "r") as weka_results:
		for i in xrange(5):
			weka_results.next()

		for i, row in enumerate(weka_results):
			if "+" not in row and row != "\n":
				conf = float(row.split()[-1])
				results[i] =  conf
	

	amines_results = [] 
	raw_rxn = DRP.models.convert_Data_to_list(reaction)
	for i, amine in enumerate(amine_list):
		best_candidate, best_conf = get_candidates(results,i*steps_per_amine, raw_rows)
		if best_candidate:
			amines_results.append( (best_conf, sim(best_candidate, raw_rxn), best_candidate) )


	if debug: print amines_results
	amines_results.sort(key=lambda x: x[0]*x[1], reverse=True)
	return amines_results

def constructRecsFromSeed(seed_pid):
  max_recommendations = 75

  #Load the default amines.
  try:
    amine_lab = DRP.models.Lab_Group.objects.filter(lab_title = "default_amines")[0]
    amines_raw = DRP.models.CompoundEntry.objects.filter(lab_group = amine_lab)
    amines_names = [c.compound for c in amines_raw]
  except Exception as e:
    raise Exception("Could not load default_amines...\n{}".format(e))

  #Get the datum and lab group
  try:
    rxn = DRP.models.Data.objects.get(id=seed_pid)
    lab = rxn.lab_group
  except Exception as e:
    raise Exception("Could not use Datum \"{}\" as seed...\n{}".format(seed_pid, e))

  #Translate abbrevs to full compound names so they play well with Paul's scripts.
  try:
    for field in reactant_fields:
      abbrev = getattr(rxn, field)
      if abbrev == "":
        continue
      compound = DRP.models.CompoundEntry.objects.filter(abbrev=abbrev)[0].compound
      setattr(rxn, field, compound)

    #Actually generate the recommendations.
  except Exception as e:
    raise Exception("Could not translate all abbrevs to compounds...\n{}".format(e))

  #Actually create the recommendations from the supplied amines and Datum.
  try:
    recommendation_list = generate_grid(rxn, amines_names)
  except Exception as e:
    raise Exception("Could not generate_grid for Datum: {}\n{}".format(seed_pid, e))

  if not recommendation_list:
    raise Exception("No recommendations generated from Datum: {}".format(seed_pid))

  #Only return a number of the possible recommendations.
  return recommendation_list[:max_recommendations]


#Generate seeds from a Datum's id.
if __name__ == "__main__":
  try:
    print constructRecsFromSeed(sys.argv[1])
    print "Seed Recommendations Built Successfully!"
  except Exception as e:
    print e
    print "\nUSAGE: python this_script.py <datum.id>"
