import json, os, sys

def process_line(line):
	try:
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
				molecules.append(abbrev_to_IUPAC[molname])
				names.append(molname)
		if oxalate:
			return ""

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
	except Exception as e:
		sys.stderr.write(e)
		return ""
	return name

if __name__ == "__main__":
	nodes = []
	abbrev_to_type = None
	abbrev_to_IUPAC = None
	analyzed = None
	
	CHEMML_DIR = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
	PREFIX_DIR = "{0}/{1}/".format(CHEMML_DIR, sys.argv[1])	
	DOT_DIR = "{0}/{1}/dots/".format(CHEMML_DIR, sys.argv[1])	
	
	restart_json_path = "{0}/restart.json".format(PREFIX_DIR)
	analyzed_json_path = "{0}/analyzed.json".format(PREFIX_DIR)
	new_nodes_txt_path = "{0}/new_nodes.txt".format(PREFIX_DIR)
	colored_json_path = "{0}/colored.json".format(PREFIX_DIR)
	edges_dot_path = "{0}/edges.dot".format(DOT_DIR)
	
	with open(new_nodes_txt_path) as node_json:
		nodes = json.load(node_json)
	with open(restart_json_path) as restart_json:
		restart = json.load(restart_json)
		abbrev_to_type = restart['abbrev_to_type']
		abbrev_to_IUPAC = restart['abbrev_to_IUPAC']
	with open(analyzed_json_path) as analyzed_json:
		analyzed = json.load(analyzed_json)
	node_set = set(nodes.keys())
	for analyzee in analyzed.keys():
		total = 0
		success = 0
		for k in analyzed[analyzee]['time']:
			total += analyzed[analyzee]['time'][k][2]
			success += analyzed[analyzee]['time'][k][0]
		if float(total) != 0 and float(success)/float(total) > 0.6:
			if analyzee in node_set:
				if nodes[analyzee] == "blue":
					nodes[analyzee] = "green"
				else:
					nodes[analyzee] = "orange"
			else:
				node_set.add(analyzee)
				nodes[analyzee] = "red"

	with open(colored_json_path,"w") as colored:
		json.dump(nodes, colored)
	with open(edges_dot_path, "a") as edges:
		for n in nodes:
			edges.write("{} [style=filled,color={}]\n".format(n, nodes[n]))
		edges.write("}")
