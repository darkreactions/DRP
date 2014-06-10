#Grab the Django settings if they aren't already set.
import sys, os, json
django_dir = os.path.dirname(os.path.realpath(__file__)).split("DRP")[0]
django_path = "{}/DRP".format(django_dir)
if django_path not in sys.path:
  sys.path.append("{}/DRP".format(django_dir))

os.environ['DJANGO_SETTINGS_MODULE'] = 'DRP.settings'


def get_cg():
	from DRP.models import CompoundEntry
	import json



	cg = dict()
	for compound in CompoundEntry.objects.all():
		if compound.calculations:
			#cg[compound.abbrev] = json.loads(compound.calculations.json_data)
			cg[compound.compound] = json.loads(compound.calculations.json_data)
		else:
			print "No calculations: {0}".format(compound.compound)

	return cg


if __name__ == "__main__":
	print "You probably meant to import this..."
