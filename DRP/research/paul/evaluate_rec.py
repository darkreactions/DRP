import json, math
from metrics import Metric

cg_props = json.load(open("restart.json"))
idxes = [0,3,6,9,12]




def evaluate_metrics(seed_set, target_set, metrics):
	seed_ckey_set = __get_ckeys(seed_set)
	target_ckey_set = __get_ckeys(target_set)

	target_ckey_set = target_ckey_set - seed_ckey_set

	possible_set = build_possible_set() - seed_ckey_set
	return [evaluate(target_ckey_set, possible_set, seed_ckey_set, m) for m in metrics ]

	
	
def evaluate(target_ckey_set, possible_set, seed_ckey_set, metric):
	scores = []


	for t in possible_set:
		scores.append((t, __compare(t, seed_ckey_set, metric)))

	from operator import itemgetter
	scores.sort(key=itemgetter(1))
	print scores
	rank = dict()
	for i in range(len(scores)):
		if scores[i][0] in target_ckey_set:
			rank[scores[i][0]] = i / float(len(possible_set))
	print rank
	return rank 


def build_possible_set():
	inorgs, orgs = split_by_type()
	print len(inorgs)*len(orgs)*len(orgs)
	gen = possible_generator(inorgs, orgs)
	ckeys = set()
	for p in gen:
		ckeys.add(__make_ckey(p, range(len(p))))
	return ckeys

def split_by_type():
	inorgs = list()
	orgs = list()
	for c in cg_props:
		if cg_props[c]["type"] == "Inorg":
			inorgs.append(c)
		elif cg_props[c]["type"] == "Org":
			orgs.append(c)
		else:
			print cg_props[c]["type"], c
	return inorgs, orgs

def possible_generator(inorgs, orgs):
	for i in range(len(inorgs)):
		for j in range(len(orgs)):
			for k in range(j+1, len(orgs)):
				yield (inorgs[i], orgs[j], orgs[k])

def __compare(target, seed_ckey_set, metric):
	best = 10000000.0
	for other in seed_ckey_set:
		best = min(best, metric.dissimilarity(target,other))
	return best


def __make_ckey(rxn, override=None):
	if override:
		l_indexes = override
	else:
		l_indexes = idxes
	ckey = set()
	for idx in l_indexes:
		ckey.add(rxn[idx])
	ckey = tuple(sorted(filter(lambda x: x != "water" and x != "", list(ckey))))
	for i in ckey:
		cg_props[i]
	return ckey


def __get_ckeys(rxns):
	ckeys = set()
	for rxn in rxns:
		try:
			ckeys.add(__make_ckey(rxn))
		except Exception as e:
			continue
	return ckeys



def __load_data(src):
	import csv
	data = []
	with open(src) as d:
		r = csv.reader(d)
		r.next()
		data = [l[1:] for l in r]
	return data

def split_data(src):
	data = __load_data(src)
	cutoff = int(len(data) * 0.7)
	seed = data[:cutoff]
	target = data[cutoff:]
	return seed, target
	
	

def avg(l):
	return sum(l)/ float(len(l))
import math
def std(l):
	a = avg(l)
	return math.sqrt(sum( [(i - a)**2 for i in l]))/float(len(l))

def calc_stats(results, name):
	def quart(s, percent):
		return s[int(percent*len(results))]
	results = results.values()
	results.sort()
	stats = dict()
	stats['name'] = name
	stats['mean'] = avg(results)
	stats['std'] = std(results)
	stats['min'] = min(results)
	stats['max'] = max(results)
	stats['q1'] = quart(results, .25) 
	stats['q2'] = quart(results, .5)
	stats['q3'] = quart(results, .75) 
	return stats
	

def main():
	print "test"
	seed, target = split_data("/home/praccugl/DRP/data.csv")

	metric_1 = Metric("tanimoto")
	metric_2 = Metric("euclidean")

	res = evaluate_metrics(seed, target, [metric_1, metric_2])

	print res[0]
	print res[1]

	print calc_stats(res[0], 'tanimoto')
	print calc_stats(res[1], 'euclidean')


if __name__ == "__main__":
	main()
