import json,math


import os, sys
full_path = os.path.dirname(os.path.realpath(__file__))+"/"
django_path = full_path[:full_path.rfind("/DRP/")]
if django_path not in sys.path:
  sys.path = [django_path] + sys.path
  os.environ['DJANGO_SETTINGS_MODULE'] = 'DRP.settings'

import DRP.model_building.model_methods as mm
from DRP.model_building import parse_rxn, load_cg
from DRP.settings import TMP_DIR
from DRP.models import ModelStats #TODO: more fully rewrite to use ModelStats objects


# Variable Setup
cg_props = load_cg.get_cg()
ml_convert = json.load(open(django_path+"/DRP/model_building/mlConvert.json"))
joint_sim = dict()
restrict_lookup = dict()
test_variables = False

def frange(start, stop, steps, int_return=False):
  # Creates a float range with n=`steps` intervals between the `start` and `stop`.
  current = start
  step = (stop-start)/steps
  while current < stop:
    value = float("{:.5f}".format(current)) # Round recommendations to a reasonable value. #TODO: STILL HERE, CASEY. Test if it works.

    if int_return: yield int(value)
    else: yield value

    current += step


def user_recommend(combinations, similarity_map, range_map):
	''' Combinations is a SET of tuples, where each tuple is a sorted
	list of the reactants tried (with any pH / temp / mass information
	stripped out. currently assuming each tuple is length 3.

	Similarity_map maps from a reactant (key) to a sorted list of
	(name, similarity) tuples where similarity = sim(key, name).

	Range map is a map from name to (min, max), where min and max are the
	largest and smallest mass to try. We can make steps in between by using
	(mix - min)/number_steps as a step size.
	'''

	explored = set(combinations) # every combination that is tried
	recommendations = [] # list of tuples. r[0] = score, r[1] = the rec
	for combination in combinations:
		recs = explore(combination, explored,similarity_map, range_map) # NOTE: Mutates explored, adds new combinations that have been tried.
		recommendations += recs
	recs_source = sorted(recommendations, key=lambda x: x[0])
	print recs_source, "source"
	recs = sorted(dissimilarity_weighting(recs_source, similarity_map), key= lambda x: x[0])
	print recs, "reweighted"
	return recs

def reweight(choice, recs, similarity_map):
	return [ dissim(choice, rec, similarity_map) for rec in recs]

def dissim(choice, rec, similarity_map, tanimoto=False):
	if not tanimoto:
		c = make_row(list(choice[1]), 0.1, 0.1, 0.1, 0.1, 1, 36, 70) #TODO: make valid
		r = make_row(list(rec[1]), 0.1, 0.1, 0.1, 0.1, 1, 36, 70) # TODO: make valid
		return (rec[0]*(1.0*euclidean_similarity(c, r)), rec[1], rec[2])
	return (rec[0]/(1.0 + tanimoto_similarity(choice, rec, similarity_map)), rec[1], rec[2])


def tanimoto_similarity(choice, rec, similarity_map):
	choice_compounds = list(choice[1])
	rec_compounds = list(rec[1])
	pairs = []
	while len(choice_compounds):
		assert(len(rec_compounds) > 0)
		choice_compound = choice_compounds.pop(0)
		max_i = -1
		max_sim = -1
		for i in range(len(rec_compounds)):
			sim = calc_similarity(choice_compound, rec_compounds[i])
			if sim > max_sim:
				max_sim = sim
				max_i = i
		pairs.append( max_sim)
		rec_compounds.pop(max_i)
	sim = sum(pairs)
	return sim



def dissimilarity_weighting(recs_unweighted, similarity_map):
	recs = []
	while len(recs_unweighted):
		choice = recs_unweighted.pop(0)
		recs.append(choice)
		recs_unweighted = reweight(choice, recs_unweighted, similarity_map)
	return recs



def do_get_evaluate_result(args):
	new_combination, range_map, combination, similarity_map = args
	score, best = evaluate_fitness(new_combination, range_map)
	score = calc_score(score, similarity_map, combination, new_combination)
	return (score, new_combination, best)

def explore(combination, explored,similarity_map, range_map):
	def args_yielder(similarity_map, combination, explored):
		for count in range(0,total_sims):
			(i1, i2, i3) = calculate_indices(count, similarity_lengths)
			new_combination = make_combination(similarity_map, combination, (i1,i2,i3))
			if new_combination in explored:
				continue
			explored.add(new_combination)
			yield (new_combination, range_map)


	similarity_lengths = [len(similarity_map[combination[0]]), len(similarity_map[combination[1]]), len(similarity_map[combination[2]])]
	total_sims = similarity_lengths[0]*similarity_lengths[1]*similarity_lengths[2]
	recs = []
	import multiprocessing
	pool = multiprocessing.Pool(processes=5)
	recs = pool.map(do_get_evaluate_result, args_yielder(similarity_map, combination, explored))
	return recs


def calculate_indices(count, list_lengths):
	indices = []
	mod = 1
	cnt = count
	for i in range(len(list_lengths)):
		mod = list_lengths[i]
		diff = cnt % mod
		indices.append(diff)
		cnt = (cnt - diff) / mod
	return indices

def make_combination(similarity_map, combination, i):
	i1,i2,i3 = i
	names = [similarity_map[combination[0]][i1][0], similarity_map[combination[1]][i2][0], similarity_map[combination[2]][i3][0]]
	return tuple(sorted(names))


def choose_center(rxns, new_combination):
	mean = [sum(r)/float(len(r)) for r in zip(*rxns)]
	dist = 100000000
	center = rxns[0]
	for rxn in range(len(rxns)):
		# calculate the distance from the mean
		d = sum( (rxns[rxn][i] - mean[i])**2 for i in range(len(rxns[rxn])))
		if d < dist: # if we're closer than any previous, set a new center
			dist = d
			center = rxns[rxn]

	return make_row(new_combination,  center[0], center[1],
		center[2], center[3], center[4], center[5], center[6])


def evaluate_fitness(new_combination, range_map, var_ranges, debug=True):
  def getDifferentReactions(tuples, num_to_get):
    def difference(rxn1, rxn2):
      weights = [0,  0,20,0,  0,20,0,  0,20,0,  0,0.5,0]
      difference = 0
      for i, weight in enumerate(weights):
        if weight>0:
          difference += abs(weight*(rxn1[i]-rxn2[i]))
      return difference

    # Only consider the tuples which have similar confidences.
    conf_threshold = 0.1
    conf = tuples[0][0]*(1.0-conf_threshold) # Assume first = max conf.
    tuples = filter(lambda tup: tup[0]>=conf, tuples)

    # Gather as many different tuples as possible.
    gathered = [tuples.pop(0)]
    diff_thresh = 0.1
    while tuples and num_to_get>0:
      tup1 = tuples.pop(0)

      if all([difference(tup1[1], tup2[1])>diff_thresh for tup2 in gathered]):
        diff = min([difference(tup1[1], tup2[1]) for tup2 in gathered])
        print "Got one:\n{} (diff: {})".format(tup1, diff)
        gathered.append(tup1)
        num_to_get -= 1
      else:
        print ".",

    return gathered

  import time, random
  from DRP.fileFunctions import createDirIfNecessary

  debug_samples = False
  return_limit=3
  search_space_max_size = 1000 #float("inf")

  # Variable and Directory Preparation.
  createDirIfNecessary(TMP_DIR)
  arff_fields, unused_indexes = mm.get_used_fields(by_current_model=True)

  # Generate different permutations of the new_combination of reactants.
  if debug:
    print "___"*10
    print "Starting row generator..."

  import sys
  #if debug:
  #  sys.stdout.write("new_combination: " + str(new_combination) + "\n")
  row_generator = generate_rows_molar(new_combination, range_map, var_ranges)
  rows = [row for row in row_generator]
  #if debug:
  #  sys.stdout.write("rows: " + str(rows) + "\n")
  #  sys.stdout.write("row len in rows: " + str(len(rows[0])) + "\n")

  # Shuffle the rows such that the search_space_max_size doesn't block combos.
  random.shuffle(rows)

  # Put the reactions in an appropriate format for handing off to WEKA by removing fields that the model doesn't know.
  cleaned = []
  for i, row in enumerate(rows[:search_space_max_size]):
    expanded = parse_rxn.parse_rxn(row, cg_props, ml_convert)
    #if debug:
    #  sys.stdout.write("expanded row len: " + str(len(expanded)) + "\n")
    cleaned.append( mm.removeUnused(expanded, unused_indexes) )
    #if debug:
    #  sys.stdout.write("cleaned row len: " + str(len(cleaned[-1])) + "\n")

  #if debug:
  #  rect = all(len(cleaned[i]) == len(cleaned[0]) for i in range(len(cleaned)))
  #  sys.stdout.write("data for arff rectangular? " + str(rect) + "\n")
  #  sys.stdout.flush()

  if debug:
    print "Search-space size: {} of {}".format(len(cleaned), len(rows))
    if debug_samples:
      print "Search-space Sample:"
      print rows[0]

  # Write all the reactions to an ARFF so that WEKA can read them.
  suffix = "_recommend"
  name = str(int(time.time()))+suffix
  if debug:
    sys.stdout.write("name of generated arff: " + repr(name) + "\n")
    #sys.stdout.write("row length of cleaned: " + str(len(cleaned[0])) + "\n")
    sys.stdout.flush()
  mm.make_arff(name, cleaned, raw_list_input=True, debug=False)

  # Run the reactions through the current WEKA model.
  if debug:
    import os
    os.system('echo "' + "about to try the path" + '"|espeak')
  current_model = ModelStats.objects.last() # TODO: CHANGE TO ModelStats.objects.filter(active=True).last() BEFORE MAKING THIS LIVE #daniel
  model_path = current_model.get_path() # used to be mm.get_current_model()

  #name = "_changed" # temporary, to see if things run with properly configured arff
  results_location = mm.make_predictions(TMP_DIR + name + ".arff", model_path, debug=debug)

  # Get the (confidence, reaction) tuples that WEKA thinks will be "successful".
  good_reactions = mm.get_good_result_tuples(results_location, rows, debug=debug)
  if debug:
    print "Good Reactions: {}".format(len(good_reactions))

  if not good_reactions:
    return [(0.0, [])]
  else:
    good_reactions.sort(key=lambda tup: tup[0])
    result = getDifferentReactions(good_reactions, return_limit)

    if debug_samples:
      print "Result Sample:"
      for row in result[:3]:
        print row

    return result




def generate_rows_molar(reactants, mass_map, var_ranges):
  def molarRange(compound, mass_map, steps):
    from DRP.compoundGuideFunctions import getMoles

    min_mass = mass_map[compound][0] if compound in mass_map else 0.1
    max_mass = mass_map[compound][1] if compound in mass_map else 0.5

    # Don't let masses get *too* extreme.
    if max_mass>8:
      max_mass = 8

    min_mols = getMoles(min_mass, compound)
    max_mols = getMoles(max_mass, compound)

    radius = 0.25

    min_mols *= (1.0-radius)
    max_mols *= (1.0+radius)
    return frange(min_mols, max_mols, steps)

  def lookupRange(field, var_ranges, var_steps, int_return=False):
    minimum = var_ranges[field][0]
    maximum = var_ranges[field][1]
    return frange(minimum, maximum, var_steps, int_return=int_return)

  def fillRow(combo, mols1, mols2, mols3, water_mols, pH, time, temp):
    from DRP.compoundGuideFunctions import getMass
    m1 = getMass(mols1, combo[0])
    m2 = getMass(mols2, combo[1])
    m3 = getMass(mols3, combo[2])
    water_mass = getMass(water_mols, "water")
    result = ["--", combo[0], m1, "g", combo[1], m2, "g", combo[2], m3, "g",
              "water", water_mass, "g", "", "","", temp, time, pH,
              "yes", "no", 4, 2,""]
    return result

  var_steps = 3
  amt_steps = 4

  # For debugging:
  import sys
  sys.stdout.write("in generate_rows_molar: reactants: " + str(reactants) + "\n")

  for pH in lookupRange("pH", var_ranges, var_steps, int_return=True):
    for time in lookupRange("time", var_ranges, var_steps, int_return=True):
      for temp in lookupRange("temp", var_ranges, var_steps, int_return=True):
        for mol1 in molarRange(reactants[0], mass_map, amt_steps):
          for mol2 in molarRange(reactants[1], mass_map, amt_steps):
            for mol3 in molarRange(reactants[2], mass_map, amt_steps):
              for water in molarRange("water", mass_map, amt_steps):
                yield fillRow(reactants, mol1, mol2, mol3, water, pH, time, temp)



def make_row(combination, m1, m2, m3, m4, pH, time, temp):
	return ["--", combination[0], m1, "g", combination[1], m2, "g", combination[2], m3, "g",  "water", m4, "g", "", "","", temp, time, pH, "yes", "no", 4, 2,""]

def calc_score(score, similarity_map, combination, new_combination):
	return score # TURNING OFF SIMILARITY TO OBSERVE EfFECT ON RECOMmeNDATIONS
	sim = []
	for i in range(len(combination)):
		for p in similarity_map[combination[i]]:
			if p[0] == new_combination[i]:
				sim.append(p[1])
				break
	from operator import mul
	return reduce(mul, sim)*score


def euclidean_similarity(reaction_one, reaction_two):
        from DRP.model_building.rxn_calculator import headers
	row_1 = parse_rxn.parse_rxn(reaction_one, cg_props, ml_convert)
	row_2 = parse_rxn.parse_rxn(reaction_two, cg_props, ml_convert)
	dist = 0
	# TODO: need to mean center and divide by std for every prop
	for i in range(len(row_1)):
		if "XXX" in headers[i]:
			continue
		if row_1[i] in ["yes", "no"]:
			if row_1[i] != row_2[i]: dist += 1
		else:
			try:
				dist += math.sqrt( (float(row_1[i]) - float(row_2[i]))**2)
			except Exception as e:
				print row_1[i], row_2[i]
				raise e
	dist = 1 / ( 1 + math.exp(- dist))
	return dist


def calc_similarity(compound_one, compound_two):
	if compound_one in joint_sim:
		if compound_two in joint_sim[compound_one]:
			return joint_sim[compound_one][compound_two]
	else:
		joint_sim[compound_one] = dict()

	if compound_two not in joint_sim:
		joint_sim[compound_two] = dict()

        if cg_props[compound_one.lower()]["type"] != cg_props[compound_one.lower()]["type"]:
		joint_sim[compound_one][compound_two] = 0.0
		joint_sim[compound_two][compound_one] = 0.0
		return 0.0

	from rdkit import DataStructs
	from rdkit.Chem.Fingerprints import FingerprintMols
	from rdkit import Chem

        mol_one = Chem.MolFromSmiles(str(cg_props[compound_one.lower()]["smiles"]))
	mol_two = Chem.MolFromSmiles(str(cg_props[compound_two.lower()]["smiles"]))
	fp_1 = FingerprintMols.FingerprintMol(mol_one)
	fp_2 = FingerprintMols.FingerprintMol(mol_two)
	similarity = DataStructs.FingerprintSimilarity(fp_1, fp_2)
	joint_sim[compound_one][compound_two] = similarity
	joint_sim[compound_two][compound_one] = similarity
	return similarity



def build_sim_list(name, cg_targets, count=5):
	if test_variables:
		count = 2
	def reweight_list(choice, others):
		return [ ( x[0], x[1], x[2]*(1.0 - calc_similarity(choice[0], x[0]))) for x in others ]
	from rdkit import DataStructs
	from rdkit.Chem.Fingerprints import FingerprintMols
	from rdkit import Chem
	mol = Chem.MolFromSmiles(str(cg_props[name]["smiles"]))
	fp = FingerprintMols.FingerprintMol(mol)
	sims = []
	for compound in cg_targets:
		if compound == name: continue
		try:
			mol2 = Chem.MolFromSmiles(str(cg_props[compound]["smiles"]))
			fp2 = FingerprintMols.FingerprintMol(mol2)
			sim = DataStructs.FingerprintSimilarity(fp, fp2)
			sims.append( (compound, sim, sim) )

		except:
			continue

	if len(sims) == 0:
		return [(name,1.0)]

	returned_list = []

	sims.sort(key=lambda x: x[1], reverse=True)
	if sims[0][1] == 0.0:
		return [(name, 1.0)]

	while len(returned_list) < count:
		choice = sims.pop(0)
		if choice[1] == 0.0:
			break
		returned_list.append(choice)
		sims = reweight_list(choice, sims)
		sims.sort(key = lambda x: x[2], reverse=True)

	if len(returned_list) > count:
		returned_list = returned_list[:count]
	return [ (x[0], x[1]) for x in returned_list]


def get_range(name):
	#TODO: calculate min, max, and then choose what upper and lower amount to take
	return (0.001*cg_props[name]["mw"], 0.01*cg_props[name]["mw"])


def build_combos(reaction_list):
	return [ tuple(sorted(filter(lambda x: x is not None and x != "water", [r[1],r[4],r[7],r[10],r[13]]))) for r in reaction_list]


def build_sim_map(compound_guide):
	inorgs,orgs,oxs = [],[],[]
	for r in compound_guide.keys():
		t = compound_guide[r]["type"]
		if t == "Org":
			orgs.append(r)
		elif t == "Inorg":
			inorgs.append(r)
		elif t == "Ox":
			oxs.append(r)
		else:
			print t

	similarity_map = {}
	for r in orgs:
		similarity_map[r] = build_sim_list(r, orgs)
	for r in inorgs:
		similarity_map[r] = build_sim_list(r, inorgs)
	for r in oxs:
		similarity_map[r] = build_sim_list(r, oxs)

def build_baseline(lab_group=None, debug=False):
	#TODO: test this stuff, rewrite to optionally accept a seed set
	# and, also, rewrite the test() method to use this.
	def average(scores):
		return sum(scores) / float(len(scores))

	def add_to_map(r_m, r):
		idxes = [1,4,7,10,13]
		for i in idxes:
			try:
				float(r[i+1])
			except:
				continue # There is no mass (ie: the mass is empty)
			if r[i] not in r_m:
				r_m[r[i]] = set()
			if float(r[i+1]) > 0:
				r_m[r[i]].add(float(r[i+1]))

        def is_numeric(elem):
          try:
            float(elem)
            return True
          except:
            return False

	from DRP.models import get_good_rxns
        from DRP.models import get_model_field_names

        current_model = ModelStats.objects.last() # TODO: CHANGE TO ModelStats.objects.filter(active=True).last() BEFORE MAKING THIS LIVE #daniel
        if not "ref" in get_model_field_names(model="Recommendation"):
            headers = ["ref"]+get_model_field_names(model="Recommendation") # Must add "ref" since "ref" is not included as a default model_field.
        else:
            headers = get_model_field_names(model="Recommendation")

	rxns = fix_abbrevs(get_good_rxns(lab_group=lab_group)[1:])

	combinations = set()
	range_map = dict()
	quality_map = dict() # A map of reactants-->reaction outcomes (ie: 1,2,3,4).

        fields_to_range = ["pH", "time", "temp"]
        indexes_to_range = [headers.index(field) for field in fields_to_range]

        var_ranges = {field:[] for field in fields_to_range}

        debug_counter = 0

        for rxn in rxns:
                r = rxn[:23]

                # Add the range variables to the data entry.
                for index, field in zip(indexes_to_range, fields_to_range):
                  if is_numeric( rxn[index] ):
                    var_ranges[field].append( float(rxn[index]) )

                reactants = filter(lambda x: x != 'water' and x != '', [ r[1], r[4], r[7], r[10], r[13]])

                if any([c not in cg_props for c in reactants]):
                  debug_counter += 1
                  continue

                q_key = tuple(sorted(reactants))
                combinations.add(q_key)
                add_to_map(range_map, r)
                if q_key not in quality_map:
                        quality_map[q_key] = []
                quality_map[q_key].append(int(r[-2]))




        # Remove the non-minimum/non-maximum masses from the range_map.
	for k in range_map:
		try:
                        radius = 0.15
                        masses = sorted([float(elem) for elem in range_map[k]])
                        minimum = masses[int(len(masses)*radius)]
                        maximum = masses[int(len(masses)*(1-radius))]
                        range_map[k] = (minimum, maximum)
		except:
			range_map[k] = (0,0)

	for q_key in quality_map:
		quality_map[q_key] = average(quality_map[q_key])


        for field, values in var_ranges.items():
          radius = 0.15
          values.sort()
          minimum = values[int(len(values)*radius)]
          maximum = values[int(len(values)*(1-radius))]
          var_ranges[field] = (minimum, maximum)


        if debug:
          print "{} of {}".format(debug_counter, len(rxns))

	return range_map, quality_map, combinations, var_ranges


def fix_abbrevs(rxns):
	abbrev_map, compound_set = get_abbrev_map()
	abbrev_map[''] = ''
	idxes = [1,4,7,10,13]
	for i in range(len(rxns)):
		for idx in idxes:
			if rxns[i][idx] in abbrev_map:
				rxns[i][idx] = abbrev_map[rxns[i][idx]]
	return rxns

def get_abbrev_map():
	from DRP.models import CompoundEntry as c

 	entries = c.objects.all()
	abbrev_map = dict()
	compound_set = set()
	for e in entries:
		abbrev_map[e.abbrev] = e.compound
		compound_set.add(e.compound)
	return abbrev_map, compound_set


def test():
	combinations = [ ("NH4VO3", "K2Cr2O7", 'pip'), ("V2O5", "H3PO3", "tmed"), ("NaVO3","SeO2","deta")]
	similarity_map = dict()
	for rxn in combinations:
		for reactant in rxn:
			if reactant not in similarity_map:
				similarity_map[reactant] = build_sim_list(reactant)
	range_map = {}
	print similarity_map, "sim"
	for name in cg_props:
		range_map[name] = get_range(name)
	range_map['water'] = (3.0, 5.0)
	print user_recommend(combinations, similarity_map, range_map)

class_map = {'Te':[u'sodium tellurite', u'hydrogen telluride', u'tellurium dioxide'], 'Se': [u'Selenium dioxide', u'selenic acid', u'selenous acid', u'Sodium selenite'], 'V':[u'sodium metavanadate', u'lithium metavanadate', u'sodium orthovanadate', u'sodium vanadium trioxide', u'potassium vanadiumtrioxide', u'Oxovanadium(2+) sulfate', u'potassium metavanadate', u'lithium vanadium trioxide', u'ammonium metavanadate', u'vanadium(V) oxide'] }




def combo_generator(seed):
        for i in seed[0]:
                for j in seed[1]:
                        for k in seed[2]:
                                yield (i,j,k)


def rank_possibilities(seed, tried):
        scorer = score_maker(tried)
	how_many_to_rate = 100 if test_variables else 500
        scores = []

        for combo in combo_generator(seed):
                scores.append( (scorer(combo), combo) )

	#Sort by the reactant tuples such that higher scores appear first.
	scores.sort(key=lambda x: x[0], reverse=True)

	results = []
	for i, score_tuple in enumerate(scores):
		if i==how_many_to_rate: break
		reweight_restrict(score_tuple[1], scores)
		scores.sort(key=lambda x: x[0], reverse=True)
		results.append(score_tuple)
        print "rank_possibilities: {}".format(len(results))
        return results

def reweight_restrict(item, scores):
	for score in scores:
		sim = score_combo(item, score[1])
		if sim > .9:
			score[0] *= 0.5


def score_maker(tried):
        def calc_score(t):
                max_score = 0
                if t in tried:
                        return 0.0
                for c in tried:
	                max_score = max(score_combo(t,c), max_score)
                if max_score == 1.0:
                        return 0.0
                return max_score

        return calc_score





def score_combo(combo_one, combo_two):
	global restrict_lookup
        def get_class(a):
                for cl in class_map:
			if a in class_map[cl]:
				return cl
                return None

	p = (combo_one, combo_two)
	if p in restrict_lookup:
		return restrict_lookup[p]

	pr = (combo_two, combo_one)
	if pr in restrict_lookup:
		return restrict_lookup[pr]

        c_one = list()
        c_two = list()
        c_three = list()
        for c in combo_one:
                if c in class_map['V']:
                        c_one.append(c)
                elif c in class_map['Se'] or c in class_map['Te']:
                        c_two.append(c)
                elif cg_props[c.lower()]['type'] == "Org":
                        c_three.append(c)
        if not (len(c_one) == len(c_two) == len(c_three) == 1):
		restrict_lookup[p] = 0.0
                return 0.0
        for c in combo_two:
                if c in class_map['V']:
                        c_one.append(c)
                elif c in class_map['Se'] or c in class_map['Te']:
                        c_two.append(c)
                elif cg_props[c.lower()]['type'] == "Org":
                        c_three.append(c)

        if not (len(c_one) == len(c_two) == len(c_three) == 2):
		restrict_lookup[p] = 0.0
                return 0.0
        score_list = []
        if c_one[0] == c_one[1]:
                score_list.append(1.0)
        else:
                score_list.append(0.5)


        if c_two[0] == c_one[0]:
                score_list.append(1.0)
        else:
                cl = [get_class(c_two[0]), get_class(c_two[1])]
                if None in cl or "V" in cl:
			print "None or V"
			restrict_lookup[pr] = 0.0
                        return 0.0
                if cl[0] == cl[1]:
                        score_list.append(0.5)
                else:
                        score_list.append(0.3)
	amine = calc_similarity(c_three[0], c_three[1])
	if amine == 1.0:
		amine = 0.0
        score_list.append(amine)
	restrict_lookup[p] = score_list[0]*score_list[1]*score_list[2]
        return score_list[0]*score_list[1]*score_list[2]



def build_diverse_org(max_results=250, debug=True):
	from random import shuffle
	orgs = [ c  for c in cg_props if cg_props[c]["type"] == "Org"]
	shuffle(orgs)

	results = [orgs[0]]
	for org in orgs[1:]:
		if len(results)==max_results: break

		max_sim = max([calc_similarity(org, r) for r in results])
		if max_sim < 0.9:
			results.append(org)

	if debug: print "Explored: {}; Retained: {}".format(len(orgs), len(results))
	return results


def recommendation_generator(use_lab_abbrevs=None, debug=False, bare_debug=False):
  import sys
  sys.stdout.write("here1\n")
  def remove_empty(rxn):
    return [field if field!="-1" else "" for field in rxn]

  from DRP.compoundGuideFunctions import translate_reactants

  # Variable Setup
  total_to_score = 1000 if not bare_debug else 50 # The number of possible combos to test.

  if debug: print "Building baseline..."
  range_map, quality_map, combinations, var_ranges = build_baseline(debug=debug)
  #if debug: sys.stdout.write("type(combinations): " + str(type(combinations)) + "\n")
  #if debug: sys.stdout.write("combinations: " + str(combinations) + "\n")

  """
  range_map = {'': (0, 0),
    u'hydrochloric acid': (0.0824, 2.0235),
    u'HIO3': (0.2149, 0.8759),
    u'R-3-aminoquinuclidine dihydrochloride': (0.1266, 0.7219),
    u"N,N'-diisopropylethylenediamine": (0.1019, 0.6559),
    ...
    }
  """

  if debug: print "Finding distinct recommendations..."
  seed = ( class_map['V'], class_map['Te'] + class_map['Se'], build_diverse_org(debug=debug) )
  #if debug: sys.stdout.write("type(seed): " + str(type(seed)) + "\n")
  #if debug: sys.stdout.write("seed: " + str(seed) + "\n")
  #if debug: sys.stdout.flush()

  if debug: print "Making abbrev_map..."
  abbrev_map, cs = get_abbrev_map()
  for i in range(len(seed)):
    for j in range(len(seed[i])):
      if seed[i][j] in abbrev_map:
        seed[i][j] = abbrev_map[seed[i][j]]

  if debug: print "Ranking possibilities..."
  #Scores: a list of (score, triple) tuples. (A score is an integer, a triple is a tuple of 3 rxn strings)
  scores = rank_possibilities(seed, combinations)
  #if debug: sys.stdout.write("scores total: " + str(scores) + "\n")

  if debug:
    import os, sys
    sys.stdout.write("here2\n")
    sys.stdout.flush()

  # A filter for redudant simplifications, using mutual information. Commenting out until fixed.
  ## if debug: print "Filtering scores (via mutual info)..."
  ## import DRP.recommendation.mutual_info as mutual_info
  ## scores = mutual_info.do_filter(scores, range_map)

  if debug: print "Scores ({} Total):".format(len(scores))
  scores = scores[:total_to_score]
  #if debug: sys.stdout.write("scores total: " + str(scores) + "\n")

  if debug: print "Rescoring..."

  for s in scores:
    try:

      reaction_tuples = evaluate_fitness(s[1], range_map, var_ranges, debug=debug)
      if not reaction_tuples:
        if debug:
          sys.stdout.write("branch 1 of recommendation generator\n")
          sys.stdout.flush()
        yield (0, None)
      else:
        for score, rxn in reaction_tuples:
          if use_lab_abbrevs:
            rxn = translate_reactants(use_lab_abbrevs, rxn, single=True)
          rxn = remove_empty(rxn)
          if debug: sys.stdout.write("branch 2 of recommendation generator\n")
          if debug: sys.stdout.write("score, s[0], score*s[0]: " + str(score) + ", " + str(s[0]) + ", " + str(score*s[0]) + "\n")
          if debug: sys.stdout.flush()

          yield (score*s[0], rxn)

    except Exception as e:
      if debug: sys.stdout.write("exception branch of recommendation generator\n")
      if debug: sys.stdout.flush()
      print "{} failed with {}: {}".format(s, e, type(e))
      yield (0, None)


def create_new_recommendations(lab_group, debug=True, bare_debug=True):
  from DRP.database_construction import store_new_Recommendation_list
  import time

  # Variable Setup
  tStart = time.clock()
  total = 0
  max_recs_per_call = 200 if not bare_debug else 5

  if debug: print "-- Creating recommendation generator..."
  scored_reactions = recommendation_generator(use_lab_abbrevs=lab_group, debug=debug, bare_debug=bare_debug)

  if debug: print "-- Storing recommmendations..."
  for (conf, rec) in scored_reactions:
    if conf>0:
      if debug: sys.stdout.write("in create_new_recommendations: conf>0\n")
      if debug: sys.stdout.flush()
      rec = map(str, rec)
      store_new_Recommendation_list(lab_group, [[conf]+rec], debug=debug)
      total += 1
      if debug: print " ... Finished #{}!".format(total)
    else:
      if debug:
        sys.stdout.write("in create_new_recommendations: conf not >0\n")
        sys.stdout.flush()

    if total>max_recs_per_call:
      break

  if debug or bare_debug:
    mins = (time.clock()-tStart)/60.0
    print "-- Rec. pipeline complete (took {:.4} minutes)!".format(mins)
    print "-- Constructed {} new recommendations".format(total)



if __name__ == "__main__":
	create_new_recommendations("Norquist Lab", debug=True, bare_debug=True)
	#print recommendation_generator()
	#print get_good_result_tuples("/home/cfalk/DevDRP/tmp/1412533505_recommend.out", [])

	#range_map, quality_map, combinations = build_baseline()
	#sim_map = build_sim_map(cg_props)
	#print sim_map
	#print user_recommend(combinations, sim_map, range_map)
