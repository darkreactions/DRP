def get_cg():
	import sys, os, json
	sys.path.append('/home/drp/web/darkreactions.haverford.edu/app/DRP')
	os.environ.setdefault('DJANGO_SETTINGS_MODULE', 'DRP.settings')
	from DRP.models import CompoundEntry



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
