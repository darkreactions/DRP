import load_data, load_cg
import numpy, math
from scipy.stats import multivariate_normal
import scipy.sparse
import time
def main():
	time_struct = dict() 
	time_struct['start'] = time.time()
	dataset = load_data.get_feature_vectors()
	time_struct['dataset'] = time.time()
	clean_dataset(dataset)
	time_struct['cleaned']  = time.time()

	time_struct['start_pdf'] = time.time()

	pdf = PDF(dataset)

	time_struct['done_pdf'] =  time.time()

	mut_inf = mutual_information(pdf, dataset)
	time_struct['done_mut'] = time.time()

	time_struct = {key: time_struct[key] - time_struct['start'] for key in time_struct}

	print mut_inf
	print time_struct

def main_triples():
	time_struct = dict()
	time_struct['start'] = time.time()
	
	dataset = load_data.get_feature_vectors_by_triple()
	load_data.collapse_triples(dataset)

	time_struct['dataset'] = time.time()

	triple_Muts = dict()
	clean_dataset_triples(dataset)

	time_struct['cleaned'] = time.time()

	for triple in dataset:
		try:
			triple_Muts[triple] = MutualInformation(dataset[triple])
		except Exception as e:
			print "Skipped {0}: {1}".format(str(triple), e)

	time_struct['done_mut'] = time.time()

	time_struct = {key: time_struct[key] - time_struct['start'] for key in time_struct }


	triple_mut_inf = dict()
	for triple in dataset:
		triple_mut_inf[triple] = triple_Muts[triple].mut_inf

	print triple_mut_inf

	max_key = triple_mut_inf.keys()[0] 
	min_key = triple_mut_inf.keys()[0] 

	for key in triple_mut_inf:
		if triple_mut_inf[key] > triple_mut_inf[max_key]:
			max_key = key
		if triple_mut_inf[key] < triple_mut_inf[min_key]:
			min_key = key

	print max_key, min_key
	print triple_mut_inf[max_key], triple_mut_inf[min_key]
	 
	print time_struct 


def get_triple_mut_from_features(dataset):
	dataset = load_data.get_feature_vectors_by_triple()
	load_data.collapse_triples(dataset)

	time_struct['dataset'] = time.time()

	triple_Muts = dict()
	clean_dataset_triples(dataset)

	time_struct['cleaned'] = time.time()

	for triple in dataset:
		try:
			triple_Muts[triple] = MutualInformation(dataset[triple])
		except Exception as e:
			print "Failed on mut {0}: {1}".format(str(triple), e)
	return triple_Muts

def find_best_row_for_delta_mut(triple_Muts, candidates):
	best_delta_mut = 0
	best_cand = -1
	best_mut = 0

	infos = []

	for i in range(len(candidates)):
		best_triple = None
		best_prob = 0.0

		for triple in triple_Muts:
			prob = triple_Muts[triple].pdf.p_feature(candidates[i])
			#print "Prob of {0} for triple {1}: {2}".format(i,triple,prob)
			if prob > best_prob:
				best_triple = triple
				best_prob = prob
		if best_triple is None:
			infos.append( {'idx':i, 'best_prob': 0, 'best_triple':None, 'delta_mut': None, 'previous_mut': None} )
			print "whaaat?"
			continue
		delta_mut = triple_Muts[best_triple].change_in_mutual(candidates[i])
		infos.append( { 'idx': i, 'best_prob':  best_prob, 'delta_mut': delta_mut, 'best_triple': best_triple, 'previous_mut': triple_Muts[best_triple].mut_inf} )
		if delta_mut > best_delta_mut:
			best_delta_mut = delta_mut
			best_cand = i
			best_mut = triple_Muts[best_triple].mut_inf

	print infos

	return best_cand, best_delta_mut, best_mut

def find_delta_mut(triple_muts, candidate):
	best_triple = None
	best_prob = 0.0
	for triple in triple_muts:
		prob = triple_muts[triple].pdf.p_feature(candidate)
		if prob > best_prob:
			best_triple = triple
			best_prob = prob
	if best_triple is None:
		return 0.0
	delta_mut = triple_muts[best_triple].change_in_mutual(candidate)
	return delta_mut

def clean_dataset_triples(dataset):
	for triple in dataset:
		clean_dataset(dataset[triple])

def clean_row(row):
	for j in range(len(row)):
		if row[j] == "no":
			row[j] = 0
		elif row[j] == "yes":
			row[j] = 1
		elif row[j] == '?' and j == 2:
			row[j] = 1
		elif row[j] == '?' and j == 4:
			row[j] = 0
		else:
			row[j] = float(row[j])

def clean_dataset(dataset):
	for i in range(len(dataset)):
		clean_row(dataset[i])


def mutual_information(pdfs, dataset):
	mut_inf = 0.0

	if not pdfs.mut_inf:
		return None 

	label = pdfs.four_mean
	not_label = pdfs.not_four_mean
	for row in dataset:
		feature = pdfs.p_feature(row)
		if feature == 0: continue

		joint = pdfs.p_joint(row, 4)
		if joint != 0:
			mut_inf += joint* math.log( joint / (label*feature), 2)

		joint = pdfs.p_joint(row, 1)
		if joint != 0:
			mut_inf += joint*math.log(joint / (not_label*feature),2)

	return mut_inf
		


class PDF:
	def __init__(self, dataset):
		self.dataset = numpy.array(dataset)
		self.labels = self.dataset[:, -1]
		self.rows = self.dataset[:, :-1]

		four_size = len(self.labels[self.labels == 4])
		if four_size <= 1 or four_size >= len(self.labels) - 1:
			raise Exception("Not enough not-4s: {0}, {1}, {2}".format(four_size, len(self.labels), self.labels.size))
		else:
			self.mut_inf = True


			self.four_mean = four_size / float(self.labels.size)
			self.not_four_mean = 1 - self.four_mean

			self.feature_mean = numpy.mean(self.rows, axis=0)
			self.feature_cov = numpy.cov(self.rows.T)
			numpy.set_printoptions(threshold=numpy.nan)


			try:
				self.feature_normal = multivariate_normal(mean=self.feature_mean, cov = self.feature_cov)
			except Exception as e:
				print dataset
				raise e

			self.cond_fours, self.cond_not_fours = make_joint_pdf(self.dataset)

	def p_outcome(self, row):
		if row[-1] == 4:
			return self.four_mean
		else:
			return self.not_four_mean

	def p_feature(self, row):
		if not self.mut_inf:
			print "womp"
			return 0.0
		row = numpy.array(row[:-1])
		return self.feature_normal.pdf(row)

	def p_joint(self, row, is_four = -1):
		if not self.mut_inf:
			return 0.0
		if is_four == -1:
			outcome = row[-1]
		else:
			outcome = is_four
		row = numpy.array(row[:-1])

		if outcome == 4:
			return self.cond_fours.pdf(row)
		else:
			return self.cond_not_fours.pdf(row)

def make_joint_pdf(dataset):
	dataset_fours = dataset[dataset[:,-1] == 4]
	dataset_fours = dataset_fours[:, :-1]

	dataset_not_fours = dataset[dataset[:, -1] != 4]
	dataset_not_fours = dataset_not_fours[:,:-1]

	mean_fours = numpy.mean(dataset_fours, axis=0)
	cov_fours = numpy.cov(dataset_fours.T)

	mean_not_fours = numpy.mean(dataset_not_fours, axis = 0)
	cov_not_fours = numpy.cov(dataset_not_fours.T)

	try:
		fours_mv = multivariate_normal(mean=mean_fours, cov = cov_fours)
	except Exception as e:
		print mean_fours, cov_fours
		print 'fours'
		raise e

	try:
		not_fours_mv = multivariate_normal(mean=mean_not_fours,cov= cov_not_fours)
	except Exception as e:
		print mean_not_fours, cov_not_fours
		print 'not fours'
		raise e
	return fours_mv, not_fours_mv


class MutualInformation:
	def __init__(self, dataset):
		self.dataset = dataset
		self.pdf = PDF(dataset)
		self.mut_inf = mutual_information(self.pdf, dataset)
		if self.mut_inf is None:
			raise Exception("No mut_inf")

	def change_in_mutual(self, row):
		new_dataset = self.dataset + [ row ] 
		new_pdfs = PDF(new_dataset)
		mut_inf = mutual_information(new_pdfs, new_dataset)
		if mut_inf is None:
			raise Exception("No new mut_inf")
		return self.mut_inf - mut_inf 

	def probability_of_row(self, row):
		return self.pdf.p_feature(row) 
		

def test_candidates():
	dataset = load_data.get_feature_vectors_by_triple()
	load_data.collapse_triples(dataset)
	triple_Muts = dict()

	clean_dataset_triples(dataset)
	for triple in dataset:
		try:
			triple_Muts[triple] = MutualInformation(dataset[triple])
		except Exception as e:
			print "failed on {0}: {1}".format(str(triple), e)

	import tmp_recs
	candidates = load_data.convert_to_feature_vectors([r[1] for r in tmp_recs.recs])
	clean_dataset(candidates)
	print  find_best_row_for_delta_mut(triple_Muts, candidates)

	

def do_filter(candidate_triples, range_map):
	import json
	mut_calc = build_mutual_calc()
	results = []

	cg = load_cg.get_cg()
	ml_convert = json.load(open("mlConvert.json"))

	for i in range(len(candidate_triples)):
		try:
			row = build_row(candidate_triples[i][1], range_map, cg, ml_convert)
			if abs(mut_calc(row)) > 0.0:
				results.append(candidate_triples[i])
			else:
				print "discarding {0}".format(candidate_triples[i])
		except Exception as e:
			print "skipping {0} due to exception: {1}".format(str(candidate_triples[i]), e)
	return results	


def build_row(triple, range_map, cg, ml_convert):
	c1 = triple[0]
	m1 = (range_map[c1][1] - range_map[c1][0])*.5 + range_map[c1][0]
	c2 = triple[1]
	m2 = (range_map[c2][1] - range_map[c2][0])*.5 + range_map[c2][0]
	c3 = triple[2]
	m1 = (range_map[c3][1] - range_map[c3][0])*.5 + range_map[c3][0]
	row = ["--", c1, m1, "g", c2, m2, "g", c3, m2, "g", "water", 6, "g", "", "", "", 120, 30, 1, "yes", "no", 4, 2, ""] 

	row = load_data.convert_to_feature_vectors([row], cg, ml_convert)[0]
	clean_row(row)
	return row 
		
		


def build_mutual_calc():
	dataset = load_data.get_feature_vectors_by_triple()
	load_data.collapse_triples(dataset)
	triple_Muts = dict()
	clean_dataset_triples(dataset)
	for triple in dataset:
		try:
			triple_Muts[triple] = MutualInformation(dataset[triple])
		except Exception as e:
			print "failed to build mutual_info struct {0}: {1}".format(str(triple), e)
	return lambda x: find_delta_mut(triple_Muts, x)	




if __name__ == "__main__":
	#main()
	#main_triples()
	test_candidates()
