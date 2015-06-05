import os, sys
full_path = os.path.dirname(os.path.realpath(__file__))+"/"
django_path = full_path[:full_path.rfind("/DRP/")]
if django_path not in sys.path:
  sys.path = [django_path] + sys.path
  os.environ['DJANGO_SETTINGS_MODULE'] = 'DRP.settings'



filename = "DRP/research/casey/raw/033115_datums.txt"
atoms = ["V", "Se"]

# Used to grab the data .
from DRP.research.casey.retrievalFunctions import get_data_from_ref_file
data = get_data_from_ref_file(filename)

from DRP.preprocessors import default_preprocessor
from DRP.model_building.rxn_calculator import headers
matrix = default_preprocessor([headers] + data)

# Model seed sets.
with open(django_path+"/DRP/research/casey/raw/ms115_spawn.txt") as f:
  ms115_set = f.read().lower().replace("\n","").split(" ")
with open(django_path+"/DRP/research/casey/raw/jho148_spawn.txt") as f:
  jho148_set = f.read().lower().replace("\n","").split(" ")
with open(django_path+"/DRP/research/casey/raw/jho213_spawn.txt") as f:
  jho213_set = f.read().lower().replace("\n","").split(" ")
with open(django_path+"/DRP/research/casey/raw/jho252_spawn.txt") as f:
  jho252_set = f.read().lower().replace("\n","").split(" ")
with open(django_path+"/DRP/research/casey/raw/033115_model.txt") as f:
  model_set = f.read().lower().replace("\n","").split(" ")
with open(django_path+"/DRP/research/casey/raw/030915_intuition.txt") as f:
  int_set = f.read().lower().replace("\n","").split(" ")


#from DRP.research.casey.knn import get_research_others, average_knn_distance
#others = get_research_others()

from DRP.models import Data
from DRP.recommendation.metrics import Tanimoto
dist = Tanimoto("test")


matrix[0] = ["XXX-Seed","XXX-Se", "XXX-Te", "XXX-Intuition", "Similarity Index"] + matrix[0]
for i, elem in enumerate(data):

  if (i%100==0):
    print "{}...".format(i)

  if elem.ref.lower() in ms115_set:
    seed = "ms115.6"
  elif elem.ref.lower() in jho148_set:
    seed = "jho148.2"
  elif elem.ref.lower() in jho213_set:
    seed = "jho213.20"
  elif elem.ref.lower() in jho252_set:
    seed = "jho252.5"
  else: seed = ""

  if seed:
    new_amine = elem.reactant_fk_3
    orig_amine = Data.objects.get(ref=seed).reactant_fk_3
    print "{}\t\t{}".format(new_amine.smiles, orig_amine.smiles)
    sim = dist.calc_similarity(new_amine, orig_amine)
  else:
    sim = ""

  #sim = average_knn_distance(elem, others, 1),

  matrix[i+1] = [seed,
                 "Se" in elem.atoms and "V" in elem.atoms,
                 "Te" in elem.atoms and "V" in elem.atoms,
                 elem.ref.lower() in int_set,
                 sim
                ] + matrix[i+1]

try:
  from DRP.fileFunctions import writeCSV
  new_csv = "DRP/research/casey/raw/full_spreadsheet_with_avg_similarities.csv"
  writeCSV(new_csv, matrix)

  print "File '{}' written!".format(new_csv)

except:
  raise Exception("Make sure you run this script from the project root!")


