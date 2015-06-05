# Prepare the Django PATH.
import os, sys
full_path = os.path.dirname(os.path.realpath(__file__))+"/"
django_path = full_path[:full_path.rfind("/DRP/")]
if django_path not in sys.path:
  sys.path = [django_path] + sys.path
  os.environ['DJANGO_SETTINGS_MODULE'] = 'DRP.settings'


import csv
from DRP.models import Data
from DRP.recommendation.metrics import Tanimoto


# Hackily load the seed and spawn columns.
with open(django_path+"/DRP/research/casey/raw/spreadsheet_seeds.txt") as f:
    seeds = f.read().lower().split("\n")[1:-1]
with open(django_path+"/DRP/research/casey/raw/spreadsheet_spawn.txt") as f:
    spawns = f.read().lower().split("\n")[1:-1]


# Prepare the similarity function.
dist = Tanimoto("test")


csv_out = "DRP/research/casey/results/similarities.csv"

with open(csv_out, "w") as f:

  # Prepare the CSV file.
  writer = csv.writer(f)
  headers = ["seed", "spawn","seed-amine","spawn-amine","seed-smiles","spawn-smiles","similarity","error"]
  writer.writerow(headers)

  for spawn, seed in zip(spawns, seeds):
    try:
      # Grab the reactions from the database.
      spawn_obj = Data.objects.get(ref__iexact=spawn)
      seed_obj = Data.objects.get(ref__iexact=seed)

      # Some of the amines are reactant 4 rather than reactant 3.
      if spawn_obj.reactant_fk_3.compound == "water":
        spawn_amine = spawn_obj.reactant_fk_4
      else:
        spawn_amine = spawn_obj.reactant_fk_3

      seed_amine = seed_obj.reactant_fk_3

      # Calculate the similarity index.
      sim = dist.calc_similarity(spawn_amine, seed_amine)

      # Prepare a CSV row.
      row = [seed, spawn, seed_amine.compound, spawn_amine.compound, seed_amine.smiles,
      spawn_amine.smiles, sim,""]

    except Exception as e:
      row = [seed, spawn, -1,-1,-1,-1,-1,"{} or {} reaction not found".format(spawn, seed)]

    writer.writerow(row)


