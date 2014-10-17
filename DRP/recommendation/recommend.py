import uuid,json,subprocess,math


import os, sys
full_path = os.path.dirname(os.path.realpath(__file__))+"/"
django_path = full_path[:full_path.rfind("/DRP/")]
if django_path not in sys.path:
  sys.path = [django_path] + sys.path
  os.environ['DJANGO_SETTINGS_MODULE'] = 'DRP.settings'

import DRP.model_building.model_methods as mm
from DRP.model_building import parse_rxn, load_cg, clean2arff
from DRP.recommendation import rebuildCDT
from DRP.settings import TMP_DIR, MODEL_DIR, BASE_DIR

cg_props = load_cg.get_cg()
ml_convert = json.load(open(django_path+"/DRP/model_building/mlConvert.json"))

joint_sim = dict()

restrict_lookup = dict()


def frange(start, stop, steps):
  current = start
  step = (stop-start)/steps
  while current < stop:
    yield current
    current += step


test_variables = True
if not test_variables:
	steps = 3
	time_range = xrange(24, 48, steps)
	pH_range = xrange(1,14, steps)
	temp_range = xrange(80, 130, steps)
else:
	steps = 3
	time_range = frange(30, 60, steps)
	pH_range = frange(1,15, steps)
	temp_range = frange(70, 100, steps)


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
	pool = multiprocessingPool(processes=5)
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


def evaluate_fitness(new_combination, range_map, debug=True):
  def removeUnused(row, unused_indexes):
    return [row[i] for i in xrange(len(row)) if i not in unused_indexes]
  def getDifferentReactions(tuples, num_to_get):
    def difference(rxn1, rxn2):
      weights = [0,  0,10,0,  0,10,0,  0,10,0,  0,0.5,0]
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
    diff_thresh = 0.01
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

  import time, csv, random

  debug_samples = False
  return_limit=3
  search_space_max_size = 5000 #float("inf")

  # Variable and Directory Preparation.
  mm.create_dir_if_necessary(TMP_DIR)
  arff_fields, unused_indexes = mm.get_used_fields()

  """
  with open("{}.csv".format(csvFilename),"w") as f:
    writer = csv.writer(f)
    writer.writerow(rebuildCDT.headers)
    for row in rows_generator:
      calcs = parse_rxn.parse_rxn(row, cg_props, ml_convert)
      writer.writerow([str(c).replace(",","c") for c in calcs])

  #clean2arff.clean(csvFilename)
  #move = "cd {};".format(django_path)
  #args = " {1} {2}".format(name, mm.get_current_model())
  #cmd = "sh {0}/DRP/model_building/make_predictions.sh".format(BASE_DIR)
  #raw_results = subprocess.check_output(move + cmd + args, shell=True)
  """



  # Generate different permutations of the new_combination of reactants.
  print "___"*10
  print "Starting row generator..."
  row_generator = generate_rows_molar(new_combination, range_map)
  rows = [row for row in row_generator]

  # Shuffle the rows such that the search_space_max_size doesn't block combos. 
  random.shuffle(rows)

  print "Calculating..."
  # Expand each row to have all the used features.
  calc_rows = [parse_rxn.parse_rxn(row, cg_props, ml_convert) for i, row in enumerate(rows)]

  # Put the reactions in an appropriate format for handing off to WEKA by removing fields that the model doesn't know.
  cleaned = [removeUnused(row, unused_indexes) for i, row in enumerate(calc_rows) if i<search_space_max_size]


  if debug:
    print "Search-space size: {} (limited to {})".format(len(rows), len(cleaned))
    if debug_samples:
      print "Search-space Sample:"
      print rows[0]

  # Write all the reactions to an ARFF so that WEKA can read them.
  suffix = "_recommend"
  name = str(int(time.time()))+suffix
  mm.make_arff(name, cleaned, raw_list_input=True, debug=False) #TODO: True=works?

  # Run the reactions through the current WEKA model.
  model_path = MODEL_DIR+mm.get_current_model()
  results_location = mm.make_predictions(TMP_DIR + name + ".arff", model_path, debug=debug)

  # Get the (confidence, reaction) tuples that WEKA thinks will be "successful".
  good_reactions = get_good_result_tuples(results_location, rows, debug=debug)
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
 
  """
  if "EMPTY" in raw_results:
    return  (0.0, get_rxn_row(0, new_combination, range_map))
  else:
    results = raw_results.split()
    print results
    result = [results[0], results[1].split(",")]
  """

  # I can reverse engineer it!
  """
  if len(result) != 2:
    raise Exception("Failed to check output: {0}".format(str(result)))
  conf, rxn_idxes = result
  """

  #rxns = [get_rxn_row(int(rxn_idx) - 1, new_combination, range_map) for rxn_idx in rxn_idxes]
  # -1 because weka isn't zero indexed.

  #rxn = choose_center(rxns, new_combination)

  #cmd = "sh /home/drp/research/chemml-research-streamlined/scripts/test_prior.sh {0}".format(name)

  #result = subprocess.check_output(cmd, shell=True)

  """
  #if float(result) < 0.2:
  #  print "prior does not like"
  #  return (-1.0, [])

  return (float(conf),rxn)
  """





def get_rxn_row(cnt, new_combination, range_map):
	ranges = range_map[new_combination[0]], range_map[new_combination[1]], range_map[new_combination[2]], range_map['water']
	mass_base = [(ranges[i][0], (ranges[i][1] - ranges[i][0])/float(steps)) for i in range(len(ranges))]
	len_list = [steps, steps, steps, steps, len(pH_range), len(time_range), len(temp_range)]
	(mass1, mass2, mass3, mass4, pH, time, temp) = calculate_indices(i, len_list)
	return [mass_base[0][0] + mass_base[0][1]*mass1, mass_base[1][0] + mass_base[1][1]*mass2, mass_base[2][0] + mass_base[2][1]*mass3, mass_base[3][1]*mass4 +mass_base[3][0], pH_range[pH], time_range[time], temp_range[temp]]


def generate_rows_molar(reactants, mass_map):
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
  for pH in frange(1,14, var_steps): 
    for time in frange(24, 48, var_steps):
      for temp in frange(80, 130, var_steps):
        for mol1 in molarRange(reactants[0], mass_map, amt_steps):
          for mol2 in molarRange(reactants[1], mass_map, amt_steps):
            for mol3 in molarRange(reactants[2], mass_map, amt_steps):
              for water in molarRange("water", mass_map, amt_steps):
                yield fillRow(reactants, mol1, mol2, mol3, water, pH, time, temp)
                

def generate_rows(new_combination, range_map):
	from operator import mul
	ranges = [range_map['water']]
	for comb in new_combination:
		if comb in range_map:
			ranges.append(range_map[comb])
		else:
			range_map[comb] = [0.001, 0.005]
			ranges.append( [0.001, 0.005] )

	# Ranges:[(0.0064, 32.995), (0.0151, 0.94), (0.1671, 2.2604), [0.001, 0.005]]
	mass_base = [(ranges[i][0], (ranges[i][1] - ranges[i][0])/float(steps)) for i in range(len(ranges))]
	len_list = [steps, steps, steps, steps, len(pH_range), len(time_range), len(temp_range)]
	#print reduce(mul, len_list)
	for i in range(reduce(mul, len_list)):
		(mass1, mass2, mass3, mass4, pH, time, temp) = calculate_indices(i, len_list)
		yield make_row(new_combination, mass_base[0][0] + mass_base[0][1]*mass1, mass_base[1][0] + mass_base[1][1]*mass2, mass_base[2][0] + mass_base[2][1]*mass3, mass_base[3][1]*mass4 +mass_base[3][0], pH_range[pH], time_range[time], temp_range[temp])


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
	row_1 = parse_rxn.parse_rxn(reaction_one, cg_props, ml_convert)
	row_2 = parse_rxn.parse_rxn(reaction_two, cg_props, ml_convert)
	dist = 0
	# TODO: need to mean center and divide by std for every prop
	for i in range(len(row_1)):
		if "XXX" in rebuildCDT.headers[i]:
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

	if cg_props[compound_one]["type"] != cg_props[compound_one]["type"]:
		joint_sim[compound_one][compound_two] = 0.0
		joint_sim[compound_two][compound_one] = 0.0
		return 0.0

	from rdkit import DataStructs
	from rdkit.Chem.Fingerprints import FingerprintMols
	from rdkit import Chem

	mol_one = Chem.MolFromSmiles(str(cg_props[compound_one]["smiles"]))
	mol_two = Chem.MolFromSmiles(str(cg_props[compound_two]["smiles"]))
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

		except Exception as e:
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
	return similarity_map

def build_baseline(lab_group=None):
	#TODO: test this stuff, rewrite to optionally accept a seed set
	# and, also, rewrite the test() method to use this.
	def quality_metric(scores):
		return sum(scores) / float(len(scores))

	def add_to_map(r_m, r):
		idxes = [1,4,7,10,13]
		for i in idxes:
			if r[i] not in r_m:
				r_m[r[i]] = set()
			try:
				float(r[i+1])
			except Exception as e:
				continue
			if float(r[i+1]) > 0:
				r_m[r[i]].add(float(r[i+1]))

	from DRP.models import get_good_rxns
	rxns = fix_abbrevs(get_good_rxns(lab_group=lab_group)[1:])

	combinations = set()
	range_map = dict()
	quality_map = dict()
	for rxn in rxns:
		r = rxn[:23]
		compoundss = filter(lambda x: x != 'water' and x != '', [ r[1], r[4], r[7], r[10], r[13]])
		if any([c not in cg_props for c in compoundss]):
			continue
		q_key = tuple(sorted(compoundss))
		combinations.add(q_key)
		add_to_map(range_map, r)
		if q_key not in quality_map:
			quality_map[q_key] = []
		quality_map[q_key].append(int(r[-2]))

	for k in range_map:
		try:
			range_map[k] = (min(range_map[k]), max(range_map[k]))
		except Exception as e:
			range_map[k] = (0,0)
	for q_key in quality_map:
		quality_map[q_key] = quality_metric(quality_map[q_key])
	return range_map, quality_map, combinations


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
                elif cg_props[c]['type'] == "Org":
                        c_three.append(c)
        if not (len(c_one) == len(c_two) == len(c_three) == 1):
		restrict_lookup[p] = 0.0
                return 0.0
        for c in combo_two:
                if c in class_map['V']:
                        c_one.append(c)
                elif c in class_map['Se'] or c in class_map['Te']:
                        c_two.append(c)
                elif cg_props[c]['type'] == "Org":
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

def recommendation_generator(use_lab_abbrevs=None, debug=False):
  def remove_empty(rxn):
    return [field if field!="-1" else "" for field in rxn]
      
  from DRP.compoundGuideFunctions import translate_reactants

  # Variable Setup
  total_to_score = 250

  if debug: print "Building baseline..."

  range_map, quality_map, combinations = build_baseline()

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

  if debug: print "Making abbrev_map..."
  abbrev_map, cs = get_abbrev_map()
  for i in range(len(seed)):
    for j in range(len(seed[i])):
      if seed[i][j] in abbrev_map:
        seed[i][j] = abbrev_map[seed[i][j]]

  if debug: print "Ranking possibilities..."
  scores = rank_possibilities(seed, combinations)

  if debug: print "Filtering scores..."
  #import DRP.research.mutual_info as mutual_info #TODO
  #scores = mutual_info.do_filter(scores, range_map)

  if debug: print "Scores ({} Total):".format(len(scores))
  scores = scores[:total_to_score]

  if debug: print "Rescoring..."
  rescored = []
  for s in scores:
    try:
      
      reaction_tuples = evaluate_fitness(s[1], range_map, debug=debug)
      if not reaction_tuples:
        yield (0, None)
      else:
        for score, rxn in reaction_tuples:
          if use_lab_abbrevs:
            rxn = translate_reactants(use_lab_abbrevs, rxn, single=True)
          rxn = remove_empty(rxn)

          yield (score*s[0], rxn)

    except Exception as e:
      print "{} failed with {}: {}".format(s, e, type(e))
      yield (0, None)


def create_new_recommendations(lab_group, debug=True):
  from DRP.database_construction import store_new_Recommendation_list
  import time
  tStart = time.clock()

  if debug: print "Creating recommendation generator..."
  scored_reactions = recommendation_generator(use_lab_abbrevs=lab_group, debug=debug)

  if debug: print "Storing recommmendations..."
  for (conf, rec) in scored_reactions:
    if conf>0:
      store_new_Recommendation_list(lab_group, [[conf]+rec], debug=debug)

  if debug: print "Recommendation construction complete (took {} minutes)!".format(time.clock()-tStart)
  


if __name__ == "__main__":
	create_new_recommendations("Norquist Lab", debug=True)
	print "COMPLETED!"
	#print recommendation_generator()
	#print get_good_result_tuples("/home/cfalk/DevDRP/tmp/1412533505_recommend.out", [])

	#range_map, quality_map, combinations = build_baseline()
	#sim_map = build_sim_map(cg_props)
	#print sim_map
	#print user_recommend(combinations, sim_map, range_map)
