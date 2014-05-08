import math, json
import load_cg
global_cg = None


def get_cg():
	if global_cg: return global_cg
	return load_cg.get_cg() 


def Metric(name):
	if name == "euclidean":
		return Euclidean(name)
	elif name == "tanimoto":
		return Tanimoto(name)
	else:
		raise Exception("Unknown metric specified, {0}".format(name))

class Euclidean:
	def __init__(self, name, models, ml_convert = None, msm = None, ):
		self.name = name
		self.cg_props = get_cg() 
		if not ml_convert:
			ml_convert = json.load(open("/home/drp/web/darkreactions.haverford.edu/app/DRP/DRP/research/mlConvert.json"))
		self.ml_convert = ml_convert 
		self.euclid_map = dict()

                if not msm:
			msm = json.load(open("mean_std_map.json"))
		self.mean_std_map, self.center_list = msm 
		import rebuildCDT
		import parse_rxn
		self.parse_rxn = parse_rxn.parse_rxn
		self.headers = rebuildCDT.headers
		self.models = models

	def make_row(self,combination):
		return ["--", combination[0], 1.0, "g", combination[1], 1.0, "g", combination[2], 1.0, "g",  "water", 5.0,  "g", "", "","", 90, 36, 1, "yes", "no", 4, 2,""]

	def apply_center(self, row):
		new_row = []
		for i in range(len(self.headers)):
			if self.center_list[i]:
				mean, std = self.mean_std_map[self.headers[i]]
				if (mean == 0 and std > 0.0001) or ( mean != 0 and abs(std / mean) > 0.0001):
					new_row.append( (float(row[i]) - float(mean)) / float(std))
				else:
					new_row.append(float(row[i]))
			else:
				new_row.append(row[i])
		return new_row

	def parse_row(self, ckey):
		if row not in self.euclid_map:
			row = self.parse_rxn(self.make_row(list(ckey)), self.cg_props, self.ml_convert)
			self.apply_center(row)
			self.euclid_map[ckey] = row
		return self.euclid_map[ckey] 

	def dissimilarity(self, A_ckey, B_ckey):
		row_1 = parse_row(A_ckey) 
		row_2 = parse_row(B_ckey) 
		return self.distance(row_1, row_2)

	def distance(self, row_1, row_2):
		dist = 0.0

		for i in range(len(row_1)):
			if "XXX" in self.headers[i]:
				continue
			if row_1[i] in ["yes", "no"] or row_2[i] in ["yes","no"]:
				if row_1[i] != row_2[i]: dist += 1.0
			else:
				try:
					dist += (float(row_1[i]) - float(row_2[i]))**2
				except Exception as e:
					print row_1, row_2
					print i
					raise e
		dist = math.sqrt(dist)
		return dist

	def sim(self, row_1, row_2):
		row_1 = self.parse_rxn(row_1, self.cg_props, self.ml_convert)
		row_2 = self.parse_rxn(row_2, self.cg_props, self.ml_convert)
		
		return self.distance(row_1, row_2)	

	def object_distance(self, obj_1, obj_2):
		row_1 = self.parse_rxn(self.models.convert_Data_to_list(obj_1), self.cg_props, self.ml_convert)
		row_2 = self.parse_rxn(self.models.convert_Data_to_list(obj_2), self.cg_props, self.ml_convert)
		return self.distance(row_1, row_2)
		





class Tanimoto:
	def __init__(self, name):
		self.compound_smiles = dict()
		self.joint_sim = dict()
		self.weird_count = 0
		from rdkit import DataStructs
		from rdkit.Chem.Fingerprints import FingerprintMols
		from rdkit import Chem
		self.ds = DataStructs
		self.FPM = FingerprintMols
		self.Chem = Chem

	def dissimilarity(self,choice, rec):
		choice_compounds = list(choice) 
		rec_compounds = list(rec)
		pairs = []
		while len(choice_compounds):
			if not len(rec_compounds) > 0:
				return 1.0
			choice_compound = choice_compounds.pop(0)
			max_i = -1
			max_sim = -1
			for i in range(len(rec_compounds)):
				sim = self.calc_similarity(choice_compound, rec_compounds[i])
				if sim > max_sim:
					max_sim = sim
					max_i = i
			pairs.append( max_sim)
			rec_compounds.pop(max_i)
		sim = sum(pairs)/len(pairs)
		return 1.0 / (1.0 + sim)


	def calc_similarity(self,compound_one, compound_two):
		if compound_one in self.joint_sim:
			if compound_two in self.joint_sim[compound_one]:
				return self.joint_sim[compound_one][compound_two]
		else:
			self.joint_sim[compound_one] = dict()
	
		if compound_two not in joint_sim:
			self.joint_sim[compound_two] = dict()
	
		if self.cg_props[compound_one]["type"] != self.cg_props[compound_one]["type"]:
			self.joint_sim[compound_one][compound_two] = 0.0
			self.joint_sim[compound_two][compound_one] = 0.0
			return 0.0
	
		
		fp_one = self.get_fp(str(self.cg_props[compound_one]["smiles"]))
		fp_two = get_fp(str(self.cg_props[compound_two]["smiles"]))
		if fp_one is None or fp_two is None:
			similarity = 0.0
			weird_count += 1
		else:
			similarity = self.ds.FingerprintSimilarity(fp_one, fp_two)
		self.joint_sim[compound_one][compound_two] = similarity
		self.joint_sim[compound_two][compound_one] = similarity
		return similarity
	
	def get_fp(smiles):
		if smiles in self.compound_smiles:
			return self.compound_smiles[smiles]
		mol = self.Chem.MolFromSmiles(smiles)
		if mol is None:
			fp = None
		else:
			fp = self.FPM.FingerprintMol(mol)
		self.compound_smiles[smiles] = fp
		return fp

