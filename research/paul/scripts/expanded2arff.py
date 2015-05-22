import sys, csv, os
from random import shuffle

def help_info():
		print "You need arguments! Usage:\n  Regular use: script.py <source_file> <new_prefix>"
		print "If you want to disable shuffle, use: script.py <source_file> <new_prefix> 1\n  For this help: script.py help"

def parse_csv(path):
	with open(path, 'r') as f:
		reader = csv.reader(f)
		row_list = [row for row in reader]
	return row_list # Does returning from with-in a "with" block cause issues?

def handle_space(row):
	# This is the most common case
	if "no" in row[1:]: #Don't care what the header i
		 r = ["no" if c!="yes" else c for c in row[1:]]
		 r.insert(0,row[0])
		 return r
	else:
		sys.stderr.write(row[0] + " has unmanaged spaces!")
		return row

def handle_comma(row): # this isn't actually necessary, is it?
	return [c.replace('"', '').replace(',','c').replace('\n','').replace('\r','') for c in row]

def not_XXX(r):
	try:
		if r[0][0:3] == "XXX":
			return False
		else:
			return True
	except Exception,e:
		sys.stderr.write("WARNING: Column header shorter than 3 characters, or empty column")
		return True

def clean_list(row_list):
	''' Does the cleaning of the list '''
	#assert: no rows have empty column headers
	#assert(all([True if r[0] != "" else False for r in row_list]))
	intermediate = map(list, zip(*row_list)) #transpose

	# Remove the reaction ID column
	intermediate = intermediate[1:]
	
	# Remove the "XXX" prefix'd columns
	intermediate = filter(not_XXX, intermediate)
	
	# Silly formatting
	intermediate = [handle_space(r) if "" in r[1:] else r for r in intermediate]
	intermediate = [handle_comma(r) if any([False if '"' in c else True for c in r]) else r for r in intermediate]
	
	# Remove purity because we only want to learn against outcome
	
	# make sure all the "yes / no" cols are handled
	final = map(list,zip(*intermediate))
	
	# save progress
	#with open(sys.argv[2] + "xp2.csv","w") as f:
	#    f.write( "\n".join([",".join(r) for r in final]))
	return final

def remove_col(row, outcome):
	if outcome:
		return row[:-2] + ["1" if row[-1] not in ['1','2','3','4',1,2,3,4] else row[-1]]
	else:
		return row[:-2] + ["1" if row[-2] not in ['1','2',1,2] else row[-2]]

def to_arff(final, special, outcome, shuffle_model=True):
	file_append = "_out" if outcome else "_pur"
	
	#Assume that script is in scripts folder.
	CHEMML_DIR = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
	file_path = "{0}/{1}/{1}".format(CHEMML_DIR, sys.argv[1])
	
	with open(file_path + file_append + ".arff", "w") as f:
		f.write("%  COMMENT \n%  NAME, DATE\n@relation " + file_path)
		for header in final[0]:
			if header in special.keys():
				if not (outcome and header == "purity") and not (not outcome and header == "outcome"):
					f.write("\n@ATTRIBUTE " + header + " " + special[header])
			else:
				f.write("\n@ATTRIBUTE " + header + " NUMERIC")
		f.write("\n\n@DATA\n")
		final = final[1:] 
		if shuffle_model:
			shuffle(final)
		f.write("\n".join([",".join(remove_col(r,outcome)) for r in final]))
	

if __name__ == "__main__":
	specials = {"outcome": "{1,2,3,4}", "slowCool": "{yes,no}", 
			"leak": "{yes,no}", "purity": "{1,2}"}
	if sys.argv[1] in ["help","h","-h","--help","--h","-help"]:
		help_info()
	elif False: # TODO: Filter out bad sys.argv[1], sys.argv[2], and make sure the files aren't already touched, and allow notouch to be overridden
		pass
	else:
		do_shuffle = True
		if len(sys.argv) == 3:
			do_shuffle = False
		#Assume that script is in scripts folder.
		CHEMML_DIR = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
		expanded_csv_path = "{0}/{1}/{1}_expanded.csv".format(CHEMML_DIR, sys.argv[1])
	
		cleaned_list = clean_list(parse_csv(expanded_csv_path))
		to_arff(cleaned_list, specials, True, do_shuffle)
		to_arff(cleaned_list, specials, False, do_shuffle)
		print "---- expanded2arff.py: Success!"

