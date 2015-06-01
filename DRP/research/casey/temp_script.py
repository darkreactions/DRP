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


matrix[0] = ["XXX-Seed","XXX-Se", "XXX-Te", "XXX-Intuition"] + matrix[0]
for i, elem in enumerate(data):
  if elem.ref.lower() in ms115_set: seed = "ms115"
  elif elem.ref.lower() in jho148_set: seed = "jho148"
  elif elem.ref.lower() in jho213_set: seed = "jho213"
  elif elem.ref.lower() in jho252_set: seed = "jho252"
  else: seed = ""

  matrix[i+1] = [seed,
                 "Se" in elem.atoms and "V" in elem.atoms,
                 "Te" in elem.atoms and "V" in elem.atoms,
                 elem.ref.lower() in int_set,
                ] + matrix[i+1]

from DRP.fileFunctions import writeCSV
new_csv = "DRP/research/casey/raw/full_spreadsheet2.csv"
writeCSV(new_csv, matrix)

#filtered = data
#filtered = [entry for entry in data if all(atom in entry.atoms for atom in atoms)]

#successes = sum([1.0 for entry in filtered if entry.outcome in {"3", "4"}])


"""
print "All from {}: {}".format(filename, len(data))
print "Filtered for {}: {}".format(atoms, len(filtered))
print "Successes: {}".format(successes)
print "Success Rate: {}".format(successes / len(filtered))
"""

print "Done!"

