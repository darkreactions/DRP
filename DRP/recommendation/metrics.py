
import sys
import os
import math
django_dir = os.path.dirname(os.path.realpath(__file__)).split("DRP")[0]
django_path = "{}/DRP".format(django_dir)
if django_path not in sys.path:
    sys.path.append("{}/DRP".format(django_dir))
    os.environ['DJANGO_SETTINGS_MODULE'] = 'DRP.settings'

from DRP.model_building import load_cg

global_cg = None


def get_cg():
    if global_cg:
        return global_cg
    return load_cg.get_cg()


class Euclidean:

    def __init__(self):
        self.cg_props = get_cg()
        self.headers = None
        self.skipped_headers = set()
        self.euclid_map = dict()

        # Caches for the means/standard deviations for each field.
        self.field_means = dict()
        self.field_stds = dict()

    def setup(self, headers, universe, debug=True):

        self.headers = headers

        field_totals = {header: 0.0 for header in headers}

        for i, header in enumerate(headers):
            for datum in universe:
                try:
                    field_totals[header] += float(datum[i])
                except:
                    self.skipped_headers.add(header)

            self.field_means[header] = field_totals[header] / len(universe)

    def apply_center(self, row):
        if not self.headers:
            raise Exception(
                "Call the `setup` method on this metric before use!")

        new_row = []
        for i, header in enumerate(self.headers):
            if header in self.skipped_headers:
                new_row.append(None)

            mean = self.field_means[header]
            std = self.field_stds[header]
            if (mean == 0 and std > 0.0001) or (mean != 0 and abs(std / mean) > 0.0001):
                new_row.append((float(row[i]) - float(mean)) / float(std))
            else:
                new_row.append(float(row[i]))

        return new_row

    def dissimilarity(self, row_1, row_2):
        """
        Find the distance these points are after "centering".
        """

        row_1 = self.apply_center(row_1)
        row_2 = self.apply_center(row_2)
        return self.distance(row_1, row_2)

    def distance(self, row_1, row_2):
        dist = 0.0

        for i, header in enumerate(self.headers):
            if header in self.skipped_headers:
                continue

            if row_1[i] in {"yes", "no", True, False}:
                if row_1[i] != row_2[i]:
                    dist += 1.0
            else:
                try:
                    dist += (float(row_1[i]) - float(row_2[i]))**2
                except Exception as e:
                    print "-- Euclidean `distance` calculation failed."
                    raise e

        dist = math.sqrt(dist)

        return dist

    def obj_distance(self, obj_1, obj_2):
        row_1 = obj_1.get_calculations_list()
        row_2 = obj_2.get_calculations_list()
        return self.distance(row_1, row_2)


class Tanimoto:

    def __init__(self, name):
        self.compound_smiles = dict()
        self.joint_sim = dict()
        self.weird_count = 0
        from rdkit import DataStructs
        from rdkit.Chem.Fingerprints import FingerprintMols
        from rdkit import Chem
        self.ds = DataStructs
        self.FPM = FingerprintMols
        self.Chem = Chem

    def dissimilarity(self, choice, rec):
        choice_compounds = list(choice)
        rec_compounds = list(rec)
        pairs = []
        while len(choice_compounds):
            if not len(rec_compounds) > 0:
                return 1.0
            choice_compound = choice_compounds.pop(0)
            max_i = -1
            max_sim = -1
            for i in range(len(rec_compounds)):
                sim = self.calc_similarity(choice_compound, rec_compounds[i])
                if sim > max_sim:
                    max_sim = sim
                    max_i = i
            pairs.append(max_sim)
            rec_compounds.pop(max_i)
        sim = sum(pairs) / len(pairs)
        return 1.0 / (1.0 + sim)

    def calc_similarity(self, compound_one, compound_two):
        if compound_one in self.joint_sim:
            if compound_two in self.joint_sim[compound_one]:
                return self.joint_sim[compound_one][compound_two]
        else:
            self.joint_sim[compound_one] = dict()

        if compound_two not in self.joint_sim:
            self.joint_sim[compound_two] = dict()

        fp_one = self.get_fp(str(compound_one.smiles))
        fp_two = self.get_fp(str(compound_two.smiles))
        if fp_one is None or fp_two is None:
            similarity = 0.0
            self.weird_count += 1
        else:
            similarity = self.ds.FingerprintSimilarity(fp_one, fp_two)
        self.joint_sim[compound_one][compound_two] = similarity
        self.joint_sim[compound_two][compound_one] = similarity
        return similarity

    def get_fp(self, smiles):
        if smiles in self.compound_smiles:
            return self.compound_smiles[smiles]
        mol = self.Chem.MolFromSmiles(smiles)
        if mol is None:
            fp = None
        else:
            fp = self.FPM.FingerprintMol(mol)
        self.compound_smiles[smiles] = fp
        return fp


def collect_universe(debug=False):
    from DRP.retrievalFunctions import get_valid_data

    universe = list(get_valid_data())

    """
  if debug:
    import random
    print "-- Using debug mode..."
    random.shuffle(universe)
    universe = universe[:100]
  """

    universe = [d.get_calculations_list() for d in universe]

    return universe


def get_default_metric(universe=None, debug=False):
    from DRP.model_building.rxn_calculator import headers

    if debug:
        print "Preparing the default metric..."

    if not universe:
        universe = collect_universe(debug=debug)

    metric = Euclidean()
    metric.setup(headers, universe)

    if debug:
        print "-- Finished!"

    return metric.distance
