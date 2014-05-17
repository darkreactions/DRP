import json, itertools, os, sys
from rdkit import DataStructs
from rdkit.Chem.Fingerprints import FingerprintMols
from rdkit.Chem import MolFromSmiles

def start():
	IUPAC_to_type = dict()
	
	CHEMML_DIR = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
	PREFIX_DIR = "{0}/{1}/".format(CHEMML_DIR, sys.argv[1])	

	restart_json_path = "{0}/restart.json".format(PREFIX_DIR)
	smiles_json_path = "{0}/smiles.json".format(PREFIX_DIR)
	similarity_json_path = "{0}/similarity.json".format(PREFIX_DIR)
	
	with open(restart_json_path) as restart:
		r = json.loads(restart.read())
		IUPAC_to_type = {r['abbrev_to_IUPAC'][k]: r['abbrev_to_type'][k] for
				k in r['abbrev_to_IUPAC']}
	
	smiles = json.loads(open(smiles_json_path).read())
	IUPAC_to_fp = dict()
	for k in IUPAC_to_type:
		try:
			IUPAC_to_fp[k] = FingerprintMols.FingerprintMol(
					MolFromSmiles(smiles[k].encode('ascii')))
		except Exception as e:
			sys.stderr.write("%s failed: %s\n" % (k, str(e)))
	
	combos = itertools.combinations_with_replacement(IUPAC_to_type.keys(),2)
	
	organics = []
	inorganics = []
	oxalates = []
	lookup_similarity = dict()
	
	for p in combos:
		try:
			one = IUPAC_to_type[p[0]]
			two = IUPAC_to_type[p[1]]
			if one == two:
				s = DataStructs.FingerprintSimilarity(
						IUPAC_to_fp[p[0]], IUPAC_to_fp[p[1]])
				#p = (p[0],p[1],s)
				lookup_similarity['_'.join(sorted([p[0],p[1]]))] = s
				if one == 'o':
					organics.append(p)
				elif one == 'i':
					inorganics.append(p)
				elif one == 'ox':
					oxalates.append(p)
				else:
					sys.stderr.write("Wut? %s" % str(p))
		except Exception as e:
			sys.stderr.write("failed: %s\n" % str(e))
	
	with open(similarity_json_path,"w") as sim:
		json.dump(lookup_similarity, sim)

if __name__ == "__main__":
	start()
