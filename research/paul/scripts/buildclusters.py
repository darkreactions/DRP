import json, random, itertools, time, sys, os

CHEMML_DIR = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
PREFIX_DIR = "{0}/{1}/".format(CHEMML_DIR, sys.argv[1])	
DOT_DIR = "{0}/{1}/dots/".format(CHEMML_DIR, sys.argv[1])
similarity_json_path = "{0}/similarity.json".format(PREFIX_DIR)

#Global variables: ###C
sims = json.load(open(similarity_json_path))
edge_list = []

def make_sim(a, b):    
	p1 = pr_str(a[0],b[0])
	p2 = pr_str(a[0],b[1])
	p3 = pr_str(a[1],b[0])
	p4 = pr_str(a[1],b[1])
	if p1 not in sims.keys():
	#print p1
		return 'hack'
	if p2 not in sims.keys():
	#print p2
		return 'hack'
	if p3 not in sims.keys():
	#print p3 
		return 'hack'
	if p4 not in sims.keys():
	#print p4
		return 'hack'
	s1 = sims[p1]
	s2 = sims[p2]
	s3 = sims[p3]
	s4 = sims[p4]
	m = max([s1,s2,s3,s4])
	if m == s1 or m == s4:
		return (s1+s4)/2
	else:
		return (s2+s3)/2
		
def pr_str(a,b):
	if a == b:
		return "same"
	return str('_'.join(sorted([a,b])))

class Cluster:
	global sims ###C 
	global edge_list ###C 
	
	def __init__(self, cluster_one=None, cluster_two=None):
		
		if cluster_one.repr_str < cluster_two.repr_str:
			self.cluster_one = cluster_one
			self.cluster_two = cluster_two
		else:
			self.cluster_one = cluster_two
			self.cluster_two = cluster_one

		if random.randint(0,1):
			self.repr_point = self.cluster_one.repr_point
			self.repr_str = self.cluster_one.repr_str
		else:
			self.repr_point = self.cluster_two.repr_point
			self.repr_str = self.cluster_two.repr_str
		self.edge_str = '_'.join(sorted([cluster_one.repr_str,cluster_two.repr_str]))
		m_sim = make_sim(self.cluster_one.repr_point[0],self.cluster_two.repr_point[0])
		nm_sim = make_sim(self.cluster_one.repr_point[1],self.cluster_two.repr_point[1])
		if m_sim == 'hack' or nm_sim == 'hack':
			self.sim = 0
		else:
			self.sim = (m_sim + nm_sim)/2
		if self.cluster_one.repr_str > self.cluster_two.repr_str:
			self.edge = self.cluster_two.repr_str + " -- " + self.cluster_one.repr_str + " [weight=" + str(self.sim) + "]"
		else:
			self.edge = self.cluster_one.repr_str + " -- " + self.cluster_two.repr_str + " [weight=" + str(self.sim) + "]"
			if self.sim >= 0.65 and self.edge not in edge_list:
				edge_list.append(self.edge)



def clean(s):
	return filter(str.isalnum, s)

class Leaf:
	def __init__(self, point):
		self.repr_point = list(set(point))
		self.repr_str = 'k' + '_'.join(sorted(map(lambda x: clean(x.encode('ascii')),point[0]) + map(lambda x: clean(x.encode('ascii')), point[1]) ))
		self.cluster_one = None
		self.cluster_two = None


def start():
	#Global Variable Setup: ###C
	global sims
	global edge_list
	global PREFIX_DIR
	global DOT_DIR
	
	try:
		out_file = sys.argv[1]
		if out_file == "help" or out_file == "h":
			print "Usage: python script.py <out_file>"
			return
	except:
		pass

	restart_json_path = "{0}/restart.json".format(PREFIX_DIR)
	smiles_json_path = "{0}/smiles.json".format(PREFIX_DIR)
	edges_dot_path = "{0}/edges.dot".format(DOT_DIR)

	organics = [] 
	inorganics = [] 
	
	IUPAC_to_smiles = None
	with open(smiles_json_path) as smiles:
		IUPAC_to_smiles = json.load(smiles)
	
	with open(restart_json_path) as info:
		data = json.load(info)
		types = data['abbrev_to_type']
		IUPAC = data['abbrev_to_IUPAC']
		for abbrev in types.keys():
			if IUPAC[abbrev] not in IUPAC_to_smiles.keys():
				sys.stderr.write("uh? %s" % IUPAC[abbrev])
				continue
			if types[abbrev] == 'o':
				organics.append(IUPAC[abbrev])
			elif types[abbrev] == 'i':
				inorganics.append(IUPAC[abbrev])
	
	sims['same'] = 1
	
	clusters = []
	metals = list(itertools.combinations_with_replacement(inorganics, 2))
	nonmetals = [(k,k) for k in organics] 
	
	def choose_and_remove(m):
		if m:
			idx = random.randrange(len(m))
			return m.pop(idx)
		return None
	nodelist = []
	for pair in metals:
		m_cluster = []
		m_copy = []
		for nonmetal in nonmetals:
			l = Leaf((pair, nonmetal))
			m_cluster.append(l)
			m_copy.append(l)
		nodelist.append(l)
		if not m_cluster:
			sys.stderr.write("Empty?!?")
			continue
		while len(m_cluster) > 1:
			c = Cluster(choose_and_remove(m_cluster), choose_and_remove(m_cluster))
			nodelist.append(c)
			m_cluster.append(c)
		for i in range(500):
			e,f = random.choice(m_copy), random.choice(m_copy)
			if e == f:
				continue
			c =Cluster(e,f)
		
		clusters.append(m_cluster[0])
	
	print "---- Clustering..."
	while len(clusters) > 1:
		c = Cluster(choose_and_remove(clusters), choose_and_remove(clusters))
		nodelist.append(c)
		clusters.append(c)
	
	
	print "---- Throwing darts..."
	for i in range(3000):
		e1 = random.choice(nodelist) 
		e2 = random.choice(nodelist)
		if e1 == e2:
			continue
		Cluster(e1,e2)
	
	with open(edges_dot_path, "w") as out:
		out.write("graph M {\nnode [label=\"\"];\n")
		out.write('\n'.join(edge_list))
	
if __name__ == "__main__":
	start()
