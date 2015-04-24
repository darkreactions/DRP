import os, sys
full_path = os.path.dirname(os.path.realpath(__file__))+"/"
django_path = full_path[:full_path.rfind("/DRP/")]
if django_path not in sys.path:
  sys.path = [django_path] + sys.path
  os.environ['DJANGO_SETTINGS_MODULE'] = 'DRP.settings'


from DRP.research.casey.retrievalFunctions import get_data_from_ref_file

filename = "DRP/research/casey/raw/030915_intuition.txt"
atoms = ["V", "Se"]

# Used to grab the data .
data = get_data_from_ref_file(filename)

filtered = data
#filtered = [entry for entry in data if all(atom in entry.atoms for atom in atoms)]

successes = sum([1.0 for entry in filtered if entry.outcome in {"3", "4"}])



print "All from {}: {}".format(filename, len(data))
print "Filtered for {}: {}".format(atoms, len(filtered))
print "Successes: {}".format(successes)
print "Success Rate: {}".format(successes / len(filtered))


