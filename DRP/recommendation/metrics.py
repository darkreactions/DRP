
import sys, os, math, json
django_dir = os.path.dirname(os.path.realpath(__file__)).split("DRP")[0]
django_path = "{}/DRP".format(django_dir)
if django_path not in sys.path:
  sys.path.append("{}/DRP".format(django_dir))
  os.environ['DJANGO_SETTINGS_MODULE'] = 'DRP.settings'

from DRP.model_building import load_cg
from DRP.settings import BASE_DIR

global_cg = None

def get_cg():
  if global_cg: return global_cg
  return load_cg.get_cg()


class Euclidean:
  def __init__(self, msm = None):
    self.cg_props = get_cg()
    self.euclid_map = dict()

    if not msm:
      msm = json.load(open(BASE_DIR+"/DRP/research/mean_std_map.json"))
    self.mean_std_map, self.center_list = msm

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

  def dissimilarity(self, row_1, row_2):
    """
    Find the distance these points are after "centering".
    """

    row_1 = self.apply_center(row_1)
    row_2 = self.apply_center(row_2)
    return self.distance(row_1, row_2)

  def distance(self, row_1, row_2):
    dist = 0.0

    for i in range(len(row_1)):
      if row_1[i] in {"yes", "no", True, False}:
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

  def obj_distance(self, obj_1, obj_2):
    row_1 = obj_1.get_calculations_list()
    row_2 = obj_2.get_calculations_list()
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

    if compound_two not in self.joint_sim:
      self.joint_sim[compound_two] = dict()

    if self.cg_props[compound_one]["type"] != self.cg_props[compound_one]["type"]:
      self.joint_sim[compound_one][compound_two] = 0.0
      self.joint_sim[compound_two][compound_one] = 0.0
      return 0.0


    fp_one = self.get_fp(str(self.cg_props[compound_one]["smiles"]))
    fp_two = self.get_fp(str(self.cg_props[compound_two]["smiles"]))
    if fp_one is None or fp_two is None:
      similarity = 0.0
      self.weird_count += 1
    else:
      similarity = self.ds.FingerprintSimilarity(fp_one, fp_two)
    self.joint_sim[compound_one][compound_two] = similarity
    self.joint_sim[compound_two][compound_one] = similarity
    return similarity

  def get_fp(self, smiles):
    if smiles in self.compound_smiles:
      return self.compound_smiles[smiles]
    mol = self.Chem.MolFromSmiles(smiles)
    if mol is None:
      fp = None
    else:
      fp = self.FPM.FingerprintMol(mol)
    self.compound_smiles[smiles] = fp
    return fp


default_metric = Euclidean().distance

