import json, csv, sys, os

CHEMML_DIR = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
PREFIX_DIR = "{0}/{1}/".format(CHEMML_DIR, sys.argv[2])	

restart_json_path = "{0}/restart.json".format(PREFIX_DIR)
smiles_json_path = "{0}/smiles.json".format(PREFIX_DIR)
new_nodes_txt_path = "{0}/new_nodes.txt".format(PREFIX_DIR)

def clean(s):
	return filter(str.isalnum, s)


def build_name(s):
	return '_'.join(sorted(map(lambda x: clean(x.encode('ascii')), s)))
abbrev_to_type = dict()
abbrev_to_IUPAC = dict()

with open(restart_json_path) as restart:
	r = json.loads(restart.read())
	abbrev_to_type = r['abbrev_to_type']
	abbrev_to_IUPAC = {k:''.join(r['abbrev_to_IUPAC'][k].split()) for 
			k in r['abbrev_to_IUPAC']}
IUPAC_to_SMILES = None
with open(smiles_json_path) as smiles_file:
	IUPAC_to_SMILES = map(lambda s: ''.join(s.split()), set(json.load(smiles_file).keys()))

reaction_count = dict()
singlemolnames = set() 
with open(sys.argv[1]) as content:
	for line in csv.reader(content):
		try:
			prev = ''
			molecules = []
			names = []
			oxalate = False
			for i in range(1,11,2):
				if line[i] == 'water':
					continue
				if line[i] == 'x' or line[i] == '-1':
					continue
				else:
					if abbrev_to_type[line[i]] == 'ox':
						oxalate = True
						break
					molname = line[i]
					mole_IUPAC = abbrev_to_IUPAC[molname]
					if mole_IUPAC not in IUPAC_to_SMILES:
						sys.stderr.write("%s not in IUPAC_TO_SMILES???" % mole_IUPAC)
						oxalate = True
					molecules.append(mole_IUPAC)
					prev = abbrev_to_IUPAC[molname]
					names.append(molname)
			if oxalate:
				continue

			if len(molecules) < 4:
				o = 0
				o_name = ''
				io = 0
				io_name = ''
				for m in names:
					if abbrev_to_type[m] == "o":
						o += 1
						o_name = abbrev_to_IUPAC[m]
					else:
						io += 1
						io_name = abbrev_to_IUPAC[m]
				if o == 1:
					molecules.append(o_name)
				if io == 1:
					molecules.append(io_name)

			molecules.sort()
			if '_'.join(molecules) in reaction_count.keys():
				if int(line[17]) == 4:
					reaction_count[build_name(molecules)][0] += 1
				reaction_count[build_name(molecules)][1] += 1
			else:
				#(mol_to_type, untyped) = add_types(molecules)
				#if untyped:
				#    continue
				#reaction_dict['0'.join(molecules)] = mol_to_type 
				if int(line[17]) == 4:
					reaction_count[build_name(molecules)] = [1,1]
				else:
					reaction_count[build_name(molecules)] = [0,1]
		except Exception as e:
			sys.stderr.write("%s: %s" % (type(e), str(e)))

for node in reaction_count:
	if reaction_count[node][0] != 0:
		reaction_count[node] = 'blue'
		#print "%s [color=red,style=filled];" % node
	else:
		reaction_count[node] = 'yellow'
		#print "%s [color=yellow,style=filled];" % node

with open(new_nodes_txt_path,'w') as nodes:
	json.dump(reaction_count,nodes)
