import os
import sys
full_path = os.path.dirname(os.path.realpath(__file__)) + "/"
django_path = full_path[:full_path.rfind("/DRP/")]
if django_path not in sys.path:
    sys.path = [django_path] + sys.path
    os.environ['DJANGO_SETTINGS_MODULE'] = 'DRP.settings'

from DRP.model_building import load_data, load_cg
import numpy
import math
from scipy.stats import multivariate_normal
import scipy.sparse
import time


def main():
    time_struct = dict()
    time_struct['start'] = time.time()
    dataset = load_data.get_feature_vectors()
    time_struct['dataset'] = time.time()
    clean_dataset(dataset)
    time_struct['cleaned'] = time.time()

    time_struct['start_pdf'] = time.time()

    pdf = PDF(dataset)

    time_struct['done_pdf'] = time.time()

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

    time_struct = {key: time_struct[key] - time_struct['start'] for key in time_struct}

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
            # print "Prob of {0} for triple {1}: {2}".format(i,triple,prob)
            if prob > best_prob:
                best_triple = triple
                best_prob = prob
        if best_triple is None:
            infos.append({'idx': i, 'best_prob': 0, 'best_triple': None, 'delta_mut': None, 'previous_mut': None})
            print "whaaat?"
            continue
        delta_mut = triple_Muts[best_triple].change_in_mutual(candidates[i])
        infos.append({'idx': i, 'best_prob': best_prob, 'delta_mut': delta_mut, 'best_triple': best_triple, 'previous_mut': triple_Muts[best_triple].mut_inf})
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
        print "WOMP"
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
            try:
                row[j] = float(row[j])
            except:
                row[j] = 0


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
        if feature == 0:
            continue

        joint = pdfs.p_joint(row, 4)
        if joint != 0:
            mut_inf += joint * math.log(joint / (label * feature), 2)

        joint = pdfs.p_joint(row, 1)
        if joint != 0:
            mut_inf += joint * math.log(joint / (not_label * feature), 2)

    return mut_inf


class PDF:

    def __init__(self, dataset):
        self.dataset = numpy.array(dataset)
        # Try-catch around line is for debugging purposes, can be removed at other times
        try:
            self.labels = self.dataset[:, -1]  # I believe these are the outcomes
        except:
            #sys.stdout.write("self.dataset: " + str(self.dataset) + "\n")
            #sys.stdout.write("self.dataset.shape: " + str(self.dataset.shape) + "\n")
            #sys.stdout.write("rectangular?: " + str(all(len(i) == len(dataset[0]) for i in dataset)) + "\n")
            # sys.stdout.flush()
            ttype, value, traceback = sys.exc_info()
            raise ttype, value, traceback
        self.rows = self.dataset[:, :-1]  # ...and this is the rest of the data

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

            # TODO: Try catch is normally here, removed and code unindented for debugging purposes. Fix when code works. Same for print statements
            # try:

            print dataset
            print "self.feature_cov shape: " + str(self.feature_cov.shape)
            print "self.feature_cov determinant: " + str(numpy.linalg.det(self.feature_cov))
            # print "self.dataset.shape: " + str(self.dataset.shape)
            # print "self.labels" + str(self.labels)
            # print "self.dataset: " + str(self.dataset)
            # print "self.rows: " + str(self.rows)
            self.feature_normal = multivariate_normal(mean=self.feature_mean, cov=self.feature_cov)

            # except Exception as e:
            # print dataset
            #raise e

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

    def p_joint(self, row, is_four=-1):
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
    dataset_fours = dataset[dataset[:, -1] == 4]
    dataset_fours = dataset_fours[:, :-1]

    dataset_not_fours = dataset[dataset[:, -1] != 4]
    dataset_not_fours = dataset_not_fours[:, :-1]

    mean_fours = numpy.mean(dataset_fours, axis=0)
    cov_fours = numpy.cov(dataset_fours.T)

    mean_not_fours = numpy.mean(dataset_not_fours, axis=0)
    cov_not_fours = numpy.cov(dataset_not_fours.T)

    try:
        fours_mv = multivariate_normal(mean=mean_fours, cov=cov_fours)
    except Exception as e:
        print mean_fours, cov_fours
        print 'fours'
        raise e

    try:
        not_fours_mv = multivariate_normal(mean=mean_not_fours, cov=cov_not_fours)
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
        new_dataset = self.dataset + [row]
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
    print find_best_row_for_delta_mut(triple_Muts, candidates)


def do_filter(candidate_triples, range_map):
    import json
    from DRP.model_building import load_cg
    ml_convert = json.load(open(django_path + "/DRP/model_building/mlConvert.json"))

    write_debug_files = False
    if write_debug_files:
        import time
        t = int(time.time())
        with open("tmp/range_map_{}.tmp".format(t), "w") as f:
            json.dump(range_map, f)
        with open("tmp/candidate_triples_{}.tmp".format(t), "w") as f:
            json.dump(candidate_triples, f)

    print "Loading CG..."
    cg = load_cg.get_cg()

    print "Building Mutual Information Calculator"
    mutual_calc = build_mutual_calc()

    # For debugging:
    import sys
    import os
    blah = "in mutual_info.do_filter"
    sys.stdout.write(blah + "\n")
    os.system('echo "' + blah + '"|espeak')
    sys.stdout.flush()

    good, bad, failed = 0, 0, 0

    results = []
    for i, triple in enumerate(candidate_triples):
        # TODO: try-catch only removed for debugging purposes. Put it back when code works.
        # try:
        row = build_row(triple[1], range_map, cg, ml_convert)
        sys.stdout.write("row: " + str(row) + "\n")
        sys.stdout.flush()
        mutual_info = abs(mutual_calc(row))
        sys.stdout.write("mutual_info: " + str(mutual_info) + "\n")
        sys.stdout.flush()
        if mutual_info > 0.0:
            results.append(triple)
            good += 1
        else:
            bad += 1
        # except Exception as e:
        #failed += 1
        # print "Skipping '{0}' due to exception: {1}".format(triple, e)

    print "# Post-MI Results: {} ({} good; {} bad; {} failed)".format(len(results), good, bad, failed)
    return results


def build_row(triple, ranges, cg, ml_convert):
    """
    Variable Examples
    triple = (u'sodium vanadium trioxide', u'selenous acid', u'1-methylpiperazine')
    range_map = {'': (0, 0),
            u'hydrochloric acid': (0.0824, 2.0235),
            u'HIO3': (0.2149, 0.8759),
            u'R-3-aminoquinuclidine dihydrochloride': (0.1266, 0.7219),
            u"N,N'-diisopropylethylenediamine": (0.1019, 0.6559),
            ...
            }
    """
    def getMeanAmount(compound, ranges):
        import random
        try:
            minimum, maximum = ranges[compound]
            return (maximum - minimum) * .5 + minimum
        except:
            general_average = 0.3
            sign = -1 if random.random() > 0.5 else 1
            return general_average + sign * random.random() * general_average

    c1, c2, c3 = triple
    m1 = getMeanAmount(c1, ranges)
    m2 = getMeanAmount(c2, ranges)
    m3 = getMeanAmount(c3, ranges)
    water = getMeanAmount("water", ranges)

    raw_row = ["--", c1, m1, "g", c2, m2, "g", c3, m2, "g", "water", water, "g", "", "", "", 120, 30, 1, "yes", "no", 4, 2, ""]

    row = load_data.convert_to_feature_vectors([raw_row], cg, ml_convert)[0][0]
    clean_row(row)
    return row


def build_mutual_calc():
    # Returns a partial function

    def rect(arr_2d):
        return all(len(i) == len(arr_2d[0]) for i in arr_2d)

    # For debugging:
    import sys
    import os
    blah = "here in DRP recommendation mutual_info build_mutual_calc"
    sys.stdout.write(blah + "\n")
    os.system('echo "' + blah + '"|espeak')
    sys.stdout.flush()

    dataset = load_data.get_feature_vectors_by_triple()
    #sys.stdout.write("all rectangular?: " + str(all(rect(item) for item in dataset.items())) + "\n")
    load_data.collapse_triples(dataset)
    #sys.stdout.write("all rectangular?: " + str(all(rect(item) for item in dataset.items())) + "\n")
    #raise(Exception("just need to take a look"))
    triple_Muts = dict()
    clean_dataset_triples(dataset)

    # For debugging:
    os.system('echo "' + "here successfully" + '"|espeak')
    tried = 0
    failed = 0
    exc = None
    exn_cnts = {}
    syst_errs = []

    for triple in dataset:
        tried = tried + 1  # for above debugging code
        try:
            triple_Mut[triple] = MutualInformation(dataset[triple])
        # More debugging code
        # except SystemError:
        #        ttype, value, traceback = sys.exc_info()
        #        raise ttype, value, traceback
        except Exception as e:
            print "failed to build mutual_info struct '{0}': {1}".format(str(triple), e)

            # for above debugging code
            failed = failed + 1
            exc = e
            e_type = (type(e), e.message[:19])
            if e_type not in exn_cnts:
                exn_cnts[e_type] = 1
            else:
                exn_cnts[e_type] += 1
            if type(e) == SystemError:
                syst_errs.append(e.message)

    # for above debugging code
    if failed != 0:
        sys.stdout.write("tried: " + str(tried) + "; failed: " + str(failed) + "\n")
        sys.stdout.write("exception types and amounts: " + str(exn_cnts) + "\n")
        sys.stdout.write("system errors: " + str(syst_errs) + "\n")
        sys.stdout.flush()
        raise(exc)

    return lambda x: find_delta_mut(triple_Muts, x)


if __name__ == "__main__":
    # main()
    # main_triples()
    # test_candidates()
    import json
    print "Loading candidate_triples and range_map from file..."
    with open("tmp/candidate_triples_1414124946.tmp", "r") as f:
        candidates = json.load(f)[:500]
    with open("tmp/range_map_1414124946.tmp", "r") as f:
        range_map = json.load(f)
    results = do_filter(candidates, range_map)
    if results:
        print "Got results!"
    print "Results: {}".format(len(results))
