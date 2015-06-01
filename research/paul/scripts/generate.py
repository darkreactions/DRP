import json, itertools, csv, sys, os

def calculate(molecule, ratio, ratio_lookup): 
	'''gets the mass from the mole ratio. Ratios sum to 1.0, 
	so we can just multiply.'''
	return ratio_lookup[molecule] * ratio # TODO: is that right? 

def generate_mole_permutations(metal_pair, nonmetals, mole_ratios, ratio_lookup):
	combos = []
	for nonmetal in nonmetals:
		if metal_pair[0] == metal_pair[1]:
			for ratios in mole_ratios: #make a different list!
				combos.append([metal_pair[0], 
						calculate(metal_pair[0], ratios[0] + ratios[1], ratio_lookup), 
						"x", "-1",
						nonmetal,
						calculate(nonmetal, ratios[2], ratio_lookup),"x","-1","water"])
		else:
			for ratios in mole_ratios:
				combos.append([metal_pair[0], 
						calculate(metal_pair[0], ratios[0], ratio_lookup), 
						metal_pair[1], 
						calculate(metal_pair[1], ratios[1], ratio_lookup), 
						nonmetal,
						calculate(nonmetal, ratios[2], ratio_lookup),"x","-1","water"])
	return combos


def generate_row(idx, parameters, parameterstats):
	result = []
	for p in range(len(parameterstats)):
		result.append(parameters[p][idx/parameterstats[p]])
		idx = idx % parameterstats[p]
	result += ["no", 2,4,""]
	return result
	



def start():
	if len(sys.argv) == 1:
		print "Usage: python script.py <prefix>"
	prefix = sys.argv[1]
	if prefix == "h" or prefix == "help":
		print "Usage: python script.py <prefix>"
	
	CHEMML_DIR = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
	PREFIX_DIR = "{0}/{1}/".format(CHEMML_DIR, sys.argv[1])	
	prefix = "{0}/{1}/{1}".format(CHEMML_DIR, sys.argv[1])	
	
	header_row = ["Reaction number","Reactant 1 name","Reactant 1 mass (g)","Reactant 2 name","Reactant 2 mass (g)","Reactant 3 name","Reactant 3 mass (g)","Reactant 4 name","Reactant 4 mass (g)","Reactant 5 name","Reactant 5 mass (g)","Temp (C)","Time (h)","Slow cool","pH (x = unknown)","Leak","purity (0 = no data in notebook, 1 = multiphase; 2 = single phase)","outcome (0 = no data in notebook, 1 = no solid; 2 = noncrystalline/brown; 3 = powder/crystallites; 4 = large single crystals)","Notes"]
	parameters = json.load(open(PREFIX_DIR+"/parameters_used.json")) #list of field names
	molecules = json.load(open(CHEMML_DIR+"/scripts/mols.json")) # dictionary of abbrev -> (mass_to_mole)
	restart = json.load(open(PREFIX_DIR+"/restart.json"))
	abbrev_to_type = restart['abbrev_to_type'] # abbrev -> type ("i","o", other crap)
	metals = list()
	nonmetals = list()

	for mol in molecules:
		if abbrev_to_type[mol] == "i":
			metals.append(mol)
		elif abbrev_to_type[mol] == "o":
			nonmetals.append(mol)
	metals = metals[:5]
	nonmetals = nonmetals[:10]
	metal_pairs = itertools.combinations_with_replacement(metals,2)

	old_list = []

	ratio_list = [(4,0,0),(0,4,0),(0,0,4),(3,1,0),(3,0,1),(1,3,0),
			(0,3,1),(1,0,3),(0,1,3),(2,1,1,),(1,2,1),(1,1,2),(2,2,0),(2,0,2)] # this should make the correct triangle
	double_ratio_list = [(4,0,0),(0,0,4),(3,0,1),(1,0,3),(2,0,2)]
	ratio_list = [
			(float(x[0])*0.19+.08, 
				float(x[1])*0.19+0.08, 
				float(x[2])*0.19+0.08) for x in ratio_list] # this converts to mole ratios

	for metal_pair in metal_pairs:
		old_list += generate_mole_permutations(metal_pair, nonmetals, ratio_list, molecules) 

	with open(prefix + "_G_expanded.csv","w") as res:
		writer = csv.writer(res)
		writer.writerow(header_row)
		pstats = [len(p) for p in parameters][::-1]
		prunning = 1
		pstats2 = []
		for pc in pstats:
			prunning *= pc
			pstats2.append(prunning)
		pcount = pstats2[-1]
		pstats = pstats2[::-1][1:] + [1]
		mcount = 0
		for row in old_list:
			for idx in range(pcount):
				writer.writerow(["%d.%d" % (mcount, idx)] + row + generate_row(idx, parameters, pstats))
			mcount += 1
	print "---- generate.py: Success!"

if __name__ == "__main__":
	start()
