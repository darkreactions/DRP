import copy, itertools, csv, sys, json, os

def swap(k):
		if k == "4":
				return 0.84
		if k == "3":
				return 0.65
		if k == "2":
				return 0.46
		if k == "1":
				return 0.27
		if k == "0":
				return 0.08

def build_point(k):
		return (swap(k[0]), swap(k[1]), swap(k[2]))

vals = ["400","040","004","310","301","130","031","103","013", "220","202","022","211","121","112"] 
r = {k: build_point(k) for k in vals }
		
def mk_dist(a,b):
		return (a[0]-b[0])**2+(a[1]-b[1])**2+(a[2]-b[2])**2

def find_nearest(r1,r2,r3):
		dist = 1000.0
		p = "400"
		for k in r:
				d = mk_dist(r[k],(r1,r2,r3))  
				if d < dist:
						dist = d
						p = k
		return p


def make_name(s):
		return "_".join(sorted(map(lambda x: clean(x.encode('ascii')), s)))
def clean(s):
		return filter(str.isalnum, s)

def mass_structure():
	ratio_dict = {p: (0,0,0,0.0) for p in vals}
	return ratio_dict

def mass_idx(m_one, m_two, nm):
	if m_two == -1:
		m_two  = 0
	r = find_nearest(float(m_one),float(m_two), float(nm))
	return r


def process_files(prefix):
	try:
		with open("%s_H_expanded.csv" % (prefix)) as exp, open("%s_H_out.out" % (prefix)) as out, open("%s_H_pur.out" % (prefix)) as pur:
				while "inst" not in out.readline():
					pur.readline()
				pur.readline()
				exp.readline()
				for exp, out_res, pur_res in itertools.izip(csv.reader(exp), out,pur):
					yield (exp, "+" not in out_res, "+" not in pur_res, out_res.split()[-1])
	except Exception as e:
		sys.stderr.write(str(e))


def start():
	pH_ranges = []
	time_ranges = []
	temp_ranges = []
	abbrev_to_IUPAC = None
	
	CHEMML_DIR = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
	PREFIX_DIR = "{0}/{1}/".format(CHEMML_DIR, sys.argv[1])	
	prefix = "{0}/{1}/{1}".format(CHEMML_DIR, sys.argv[1])	    

	restart_json_path = "{0}/restart.json".format(PREFIX_DIR)
	parameters_json_path = "{0}/parameters_used.json".format(PREFIX_DIR)
	analyzed_json_path = "{0}/analyzed.json".format(PREFIX_DIR)
	hierarchy_json_path = "{0}/hierarchy.json".format(PREFIX_DIR)
	
	with open(restart_json_path) as restart:
		abbrev_to_IUPAC = json.load(restart)['abbrev_to_IUPAC']
	with open(parameters_json_path) as param:
		params = json.load(param)
		for t in params[1]:
			temp_ranges.append(t)
		for time in params[2]:
			time_ranges.append(time)
		for pH in params[4]:
			pH_ranges.append(pH)
	base = {'temp':{str(k):(0,0,0,0.0) for k in temp_ranges}, 'time':{str(k):(0,0,0,0.0) for k in time_ranges},
			'pH': {str(k):mass_structure() for k in pH_ranges}}
	parsed_data = dict()
	seen_combo = set()
	metal_pair_to_names = dict()
	metal_pair_seen = set()
	for (experiment, out_res, pur_res, conf) in process_files(prefix):
		inorg_one = experiment[1]
		inorg_moles_one = experiment[3]
		inorg_two = experiment[4]
		inorg_moles_two = experiment[6]
		org_name = experiment[10] 
		org_moles = experiment[12]
		temp = experiment[19].split(".")[0]
		time = experiment[20].split(".")[0]
		pH = experiment[22].split(".")[0]

		name = None
		inorg_one_IUPAC = None
		inorg_two_IUPAC = None
		org_IUPAC = None
		try:
			if inorg_two == "-1":
				inorg_two = inorg_one
			inorg_one_IUPAC = abbrev_to_IUPAC[inorg_one]
			inorg_two_IUPAC = abbrev_to_IUPAC[inorg_two]
			org_IUPAC = abbrev_to_IUPAC[org_name]
			name = make_name([inorg_one_IUPAC, inorg_two_IUPAC, org_IUPAC, org_IUPAC])
		except Exception as e:
			sys.stderr.write("Skipping %s" % str(e))
			continue

		if name not in seen_combo:
			seen_combo.add(name)
			parsed_data[name] = copy.deepcopy(base)
			metal_name = make_name([inorg_one_IUPAC, inorg_two_IUPAC])
			if metal_name not in metal_pair_seen:
				metal_pair_seen.add(metal_name)
				metal_pair_to_names[metal_name] = set()
			metal_pair_to_names[metal_name].add((name, inorg_one, inorg_two, org_name))


		mass_number = mass_idx(inorg_moles_one, inorg_moles_two, org_moles)
		
		o_temp,p_temp,t_temp, c_temp = parsed_data[name]['temp'][str(temp)]
		o_time,p_time,t_time, c_time = parsed_data[name]['time'][str(time)]
		o_pH, p_pH, t_pH, c_pH = parsed_data[name]['pH'][str(pH)][mass_number]
		if out_res:
			o_temp += 1
			o_time += 1
			o_pH += 1
		if pur_res:
			p_time += 1
			p_temp += 1
			p_pH += 1
		t_time += 1
		t_temp += 1
		t_pH += 1
		c_time += float(conf)
		c_temp += float(conf)
		c_pH += float(conf)
		parsed_data[name]['temp'][str(temp)] = (o_temp, p_temp, t_temp, c_temp)
		parsed_data[name]['time'][str(time)] = (o_time, p_time, t_time, c_time)
		parsed_data[name]['pH'][str(pH)][mass_number] = (o_pH, p_pH, t_pH, c_pH)

	with open(analyzed_json_path,'w') as dump:
		json.dump(parsed_data, dump)
	with open(hierarchy_json_path,'w') as dump:
		json.dump({k:list(metal_pair_to_names[k]) for k in metal_pair_to_names}, dump)
	print "---- analyzed.py: Success!"


if __name__ == "__main__":
	start()
