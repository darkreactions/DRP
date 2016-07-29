from copy import deepcopy

from sklearn.metrics.cluster import adjusted_mutual_info_score

from DRP.plugins.rxndescriptors.rxnhash import calculate
from DRP.recommender.AbstractRecommender import AbstractRecommender
from DRP.recommender.ReactionGenerator import ReactionGenerator
from DRP.recommender.ReactionSieve import ReactionSieve
from DRP.models import PerformedReaction, BoolRxnDescriptorValue, CatRxnDescriptorValue, CompoundQuantity, Reaction, Compound, RecommendedReaction, LabGroup

headers_to_use = ['Pu_mols_DRP_0.02', 'Tb_mols_DRP_0.02', 'Lu_mols_DRP_0.02', 'Ru_mols_DRP_0.02', 'Ga_mols_DRP_0.02', 'Au_mols_DRP_0.02', 'B_mols_DRP_0.02', 'Li_mols_DRP_0.02', 'Ca_mols_DRP_0.02', 'Hf_mols_DRP_0.02', 'Pd_mols_DRP_0.02', 'In_mols_DRP_0.02', 'Org_mw_drp/rdkit_0_Range_DRP_DRP', 'Pr_mols_DRP_0.02', 'Ho_mols_DRP_0.02', 'Na_mols_DRP_0.02', 'Cu_mols_DRP_0.02', 'reaction_pH_manual_0', 'Po_mols_DRP_0.02', 'Org_mw_drp/rdkit_0_gmean_molarity_DRP_DRP', 'Zn_mols_DRP_0.02', 'Tm_mols_DRP_0.02', 'Inorg_mw_drp/rdkit_0_Range_DRP_DRP', 'H_mols_DRP_0.02', 'reaction_temperature_manual_0', 'At_mols_DRP_0.02', 'Te_mols_DRP_0.02', 'Rh_mols_DRP_0.02', 'Yb_mols_DRP_0.02', 'Inorg_amount_count_DRP_0.02', 'Pt_mols_DRP_0.02', 'Pm_mols_DRP_0.02', 'Sb_mols_DRP_0.02', 'Al_mols_DRP_0.02', 'F_mols_DRP_0.02', 'S_mols_DRP_0.02', 'Mn_mols_DRP_0.02', 'Cr_mols_DRP_0.02', 'Nd_mols_DRP_0.02', 'Co_mols_DRP_0.02', 'Fe_mols_DRP_0.02', 'C_mols_DRP_0.02', 'Os_mols_DRP_0.02', 'Ti_mols_DRP_0.02', 'Ir_mols_DRP_0.02', 'Si_mols_DRP_0.02', 'Re_mols_DRP_0.02', 'Ar_mols_DRP_0.02', 'Solv_amount_count_DRP_0.02', 'W_mols_DRP_0.02', 'Dy_mols_DRP_0.02', 'Y_mols_DRP_0.02', 'Cd_mols_DRP_0.02', 'Ce_mols_DRP_0.02', 'As_mols_DRP_0.02', 'Ag_mols_DRP_0.02', 'La_mols_DRP_0.02', 'U_mols_DRP_0.02', 'Np_mols_DRP_0.02', 'Eu_mols_DRP_0.02', 'Er_mols_DRP_0.02', 'K_mols_DRP_0.02', 'Pa_mols_DRP_0.02', 'Ni_mols_DRP_0.02', 'pH_amount_molarity_DRP_0.02', 'Ox_amount_molarity_DRP_0.02', 'N_mols_DRP_0.02', 'Am_mols_DRP_0.02', 'Tl_mols_DRP_0.02', 'Kr_mols_DRP_0.02', 'Inorg_mw_drp/rdkit_0_Max_DRP_DRP', 'Ra_mols_DRP_0.02', 'Org_amount_count_DRP_0.02', 'Se_mols_DRP_0.02', 'Be_mols_DRP_0.02', 'Th_mols_DRP_0.02', 'Inorg_mw_drp/rdkit_0_gmean_count_DRP_DRP', 'Mg_mols_DRP_0.02', 'Rn_mols_DRP_0.02', 'Tc_mols_DRP_0.02', 'Cl_mols_DRP_0.02', 'Sn_mols_DRP_0.02', 'Inorg_mw_drp/rdkit_0_gmean_molarity_DRP_DRP', 'Sc_mols_DRP_0.02', 'Ac_mols_DRP_0.02', 'Gd_mols_DRP_0.02', 'Xe_mols_DRP_0.02', 'Hg_mols_DRP_0.02', 'Ge_mols_DRP_0.02', 'Mo_mols_DRP_0.02', 'Cs_mols_DRP_0.02', 'Sr_mols_DRP_0.02', 'O_mols_DRP_0.02', 'Org_mw_drp/rdkit_0_gmean_count_DRP_DRP', 'Nb_mols_DRP_0.02', 'I_mols_DRP_0.02', 'Ta_mols_DRP_0.02', 'Fr_mols_DRP_0.02', 'Rb_mols_DRP_0.02', 'Ox_amount_count_DRP_0.02', 'pH_amount_count_DRP_0.02', 'Zr_mols_DRP_0.02', 'He_mols_DRP_0.02', 'Ne_mols_DRP_0.02', 'reaction_time_manual_0', 'Pb_mols_DRP_0.02', 'Br_mols_DRP_0.02', 'Inorg_amount_molarity_DRP_0.02', 'V_mols_DRP_0.02', 'Org_mw_drp/rdkit_0_Max_DRP_DRP', 'Bi_mols_DRP_0.02', 'Sm_mols_DRP_0.02', 'Solv_amount_molarity_DRP_0.02', 'Ba_mols_DRP_0.02', 'Org_amount_molarity_DRP_0.02', 'P_mols_DRP_0.02']
headers_to_use = ['Pu_mols_DRP_0.02', 'Tb_mols_DRP_0.02', 'Lu_mols_DRP_0.02', 'Ru_mols_DRP_0.02', 'Ga_mols_DRP_0.02', 'Au_mols_DRP_0.02', 'B_mols_DRP_0.02', 'Li_mols_DRP_0.02', 'Ca_mols_DRP_0.02', 'Hf_mols_DRP_0.02', 'Pd_mols_DRP_0.02', 'In_mols_DRP_0.02', 'Org_mw_drp/rdkit_0_Range_DRP_DRP', 'Pr_mols_DRP_0.02', 'Ho_mols_DRP_0.02', 'Na_mols_DRP_0.02', 'Cu_mols_DRP_0.02', 'reaction_pH_manual_0', 'Po_mols_DRP_0.02', 'Org_mw_drp/rdkit_0_gmean_molarity_DRP_DRP', 'Zn_mols_DRP_0.02', 'Tm_mols_DRP_0.02', 'Inorg_mw_drp/rdkit_0_Range_DRP_DRP', 'H_mols_DRP_0.02', 'reaction_temperature_manual_0', 'At_mols_DRP_0.02', 'Te_mols_DRP_0.02', 'Rh_mols_DRP_0.02', 'Yb_mols_DRP_0.02', 'Inorg_amount_count_DRP_0.02', 'Pt_mols_DRP_0.02', 'Pm_mols_DRP_0.02', 'Sb_mols_DRP_0.02', 'Al_mols_DRP_0.02', 'F_mols_DRP_0.02', 'S_mols_DRP_0.02', 'Mn_mols_DRP_0.02', 'Cr_mols_DRP_0.02', 'Nd_mols_DRP_0.02', 'Co_mols_DRP_0.02', 'Fe_mols_DRP_0.02', 'C_mols_DRP_0.02', 'Os_mols_DRP_0.02', 'Ti_mols_DRP_0.02', 'Ir_mols_DRP_0.02', 'Si_mols_DRP_0.02', 'Re_mols_DRP_0.02', 'Ar_mols_DRP_0.02', 'Solv_amount_count_DRP_0.02', 'W_mols_DRP_0.02', 'Dy_mols_DRP_0.02', 'Y_mols_DRP_0.02', 'Cd_mols_DRP_0.02', 'Ce_mols_DRP_0.02', 'As_mols_DRP_0.02', 'Ag_mols_DRP_0.02', 'La_mols_DRP_0.02', 'U_mols_DRP_0.02', 'Np_mols_DRP_0.02', 'Eu_mols_DRP_0.02', 'Er_mols_DRP_0.02', 'K_mols_DRP_0.02', 'Pa_mols_DRP_0.02', 'Ni_mols_DRP_0.02', 'pH_amount_molarity_DRP_0.02', 'Ox_amount_molarity_DRP_0.02', 'N_mols_DRP_0.02', 'Am_mols_DRP_0.02', 'Tl_mols_DRP_0.02', 'Kr_mols_DRP_0.02', 'Inorg_mw_drp/rdkit_0_Max_DRP_DRP', 'Ra_mols_DRP_0.02', 'Org_amount_count_DRP_0.02', 'Se_mols_DRP_0.02', 'Be_mols_DRP_0.02', 'Th_mols_DRP_0.02', 'Inorg_mw_drp/rdkit_0_gmean_count_DRP_DRP', 'Mg_mols_DRP_0.02', 'Rn_mols_DRP_0.02', 'Tc_mols_DRP_0.02', 'Cl_mols_DRP_0.02', 'Sn_mols_DRP_0.02', 'Inorg_mw_drp/rdkit_0_gmean_molarity_DRP_DRP', 'Sc_mols_DRP_0.02', 'Ac_mols_DRP_0.02', 'Gd_mols_DRP_0.02', 'Xe_mols_DRP_0.02', 'Hg_mols_DRP_0.02', 'Ge_mols_DRP_0.02', 'Mo_mols_DRP_0.02', 'Cs_mols_DRP_0.02', 'Sr_mols_DRP_0.02', 'O_mols_DRP_0.02', 'Org_mw_drp/rdkit_0_gmean_count_DRP_DRP', 'Nb_mols_DRP_0.02', 'I_mols_DRP_0.02', 'Ta_mols_DRP_0.02', 'Fr_mols_DRP_0.02', 'Rb_mols_DRP_0.02', 'Ox_amount_count_DRP_0.02', 'pH_amount_count_DRP_0.02', 'Zr_mols_DRP_0.02', 'He_mols_DRP_0.02', 'Ne_mols_DRP_0.02', 'reaction_time_manual_0', 'Pb_mols_DRP_0.02', 'Br_mols_DRP_0.02', 'Inorg_amount_molarity_DRP_0.02', 'V_mols_DRP_0.02', 'Org_mw_drp/rdkit_0_Max_DRP_DRP', 'Bi_mols_DRP_0.02', 'Sm_mols_DRP_0.02', 'Solv_amount_molarity_DRP_0.02', 'Ba_mols_DRP_0.02', 'Org_amount_molarity_DRP_0.02', 'P_mols_DRP_0.02']
headers_to_use = ['reaction_temperature_manual_0', 'reaction_time_manual_0']
# TODO Get rid of this ASAP ^ but meet Sorelle's deadline

class Recommender(AbstractRecommender):
    def __init__(self, model_container, grid_params, desired_desc_dict, seed = None): # Reaction.objects.get(notes = 'Placeholder reaction')):
        """
        High level class for making recommendation. Takes a specific model, parameters given from a gridparameter object, a dictionary of desired descriptors and an optional seed.

        All Descriptor Dicts use (Descriptor, Value List) pairs, where the
        Value List is a list of any values this Descriptor should have.
        """

        print grid_params
        self.generator = ReactionGenerator(grid_params)

        # self.plausible_reactions]
        self.model_container = model_container
        self.desired_desc_dict = desired_desc_dict

        # TODO remove
        self.desc_dict = grid_params

        self.hash_MI = {}

        self.lshash = None
        self.headers = None

        # Might remove?
        self.seed = seed


    def build_hash_outcome_dict(self):
        """Build a dictionary of the form 'hash' : [outcome1, outcome2, outcome3]."""
        outcome_dict = {}
        count = 0
        for rxn in PerformedReaction.objects.all():
            count +=1
            try:
                # TODO TEST TO MAKE SURE HASH CATEGORICAL DESCRIPTOR PK DOesN'T CHANGE
                rxn_hash = CatRxnDescriptorValue.objects.get(descriptor__pk=6793,reaction=rxn).value
            except CatRxnDescriptorValue.DoesNotExist:
                calculate(rxn)
                try:
                    rxn_hash = CatRxnDescriptorValue.objects.get(descriptor__pk=6793,reaction=rxn).value
                except CatRxnDescriptorValue.DoesNotExist: # missing or incorrect reactant set up
                    continue

            try:
                reaction_outcome = BoolRxnDescriptorValue.objects.get(descriptor__pk=2,reaction=rxn)
            except BoolRxnDescriptorValue.DoesNotExist:
                continue

            if reaction_outcome:
                reaction_outcome = 1
            else:
                reaction_outcome = 0

            if rxn_hash in outcome_dict:
                outcome_dict[rxn_hash].append(reaction_outcome)
            else:
                outcome_dict[rxn_hash] = [reaction_outcome]

        to_remove = []

        for i in outcome_dict:
            if 1 not in outcome_dict[i] or 0 not in outcome_dict[i]:
                to_remove.append(i)
            else:
                pass

        for i in to_remove:
            del outcome_dict[i]

        return outcome_dict

    def get_mutual_info(self, outcome_dict):
        """Gets mutual information given hash strings and outcomes using sklearn function."""
        hashes, outcomes = zip(*list(outcome_dict.items()))
        return adjusted_mutual_info_score(hashes, outcomes)

    def hash_perform_generator(self, compound_dict, anash):
        """Determine the minimum change in information given a compound_dict and a hash to add."""
        deepdictlow = deepcopy(compound_dict)
        deepdicthigh = deepcopy(compound_dict)
        if anash in deepdictlow:
            deepdictlow[anash].append(0)
            deepdicthigh[anash].append(1)
        else:
            deepdictlow[anash] = [0]
            deepdicthigh[anash] = [1]
        return (anash, min(self.get_mutual_info(deepdictlow), self.get_mutual_info(deepdicthigh)))

    def get_best_delta_MI(self, k=100):
        """Get the compounds that lead to the k hashes with the highest guaranteed mutual information change"""
        base_dict = self.build_hash_outcome_dict()

        triples_dict = self.generator.get_all_triple_hashes()

        hash_min_delta_MI = []
        for rxn_hash in triples_dict:
            result = self.hash_perform_generator(base_dict, rxn_hash)
            hash_min_delta_MI.append(result)
            self.hash_MI[rxn_hash] = result[1]

        hash_min_delta_MI.sort(key=(lambda x: x[1]), reverse=True)

        greatest_delta_MI_compounds = []
        for a_hash in hash_min_delta_MI[:k]:
            greatest_delta_MI_compounds.append(triples_dict[a_hash[0]])

        return greatest_delta_MI_compounds

    def get_lshash(self):
        """Index all existing reactions based on specified headers into an lshash."""
        from lshash import LSHash
        headers = headers_to_use
        count = 0

        lsh = LSHash(1, len(headers))
        for i in PerformedReaction.objects.all().rows(True):
            to__index = []
            for header in headers:
                try:
                    to__index.append(i[header])
                except KeyError:
                    continue
            if len(to__index) == len(headers):
                lsh.index(to__index)
                count += 1
            else:
                pass

        print 'count', count
        self.lshash = lsh
        self.headers = headers

    def normalize_dist_dict(self, dist_dict):
        """Normalize all distances on a 0 to 1 scale"""
        min_knn_dist = min(dist_dict.values())
        max_knn_dist = max(dist_dict.values())
        split = max_knn_dist - min_knn_dist
        if split != 0:
            for key in dist_dict:
                dist_dict[key] = (dist_dict[key] - min_knn_dist) / split
        else:
            for key in dist_dict:
                dist_dict[key] = 1

    def get_distance_dict(self, number_of_neighbors=20, normalize=True):
        """Get the average distance of each plausible reaction to the number_of_neighbors nearest neighbors."""
        dist_dict = {}

        if self.lshash is None:
            self.get_lshash()

        for rxn, rxn_row in zip(self.plausible_reactions, self.plausible_reactions.rows(True)):
            point = []
            for header in self.headers:
               point.append(rxn_row[header])
            # Point : Average distance from 20 nearest neighbors
            dist_dict[rxn] = sum([query_result[1] for query_result in self.lshash.query(point, num_results=number_of_neighbors)])/number_of_neighbors

        if normalize:
            self.normalize_dist_dict(dist_dict)

        return dist_dict


    def score_candidates(self):
        """Determine the score for all plausible reactions."""

        score_dict = {}
        dist_dict = self.get_distance_dict()

        for rxn in self.plausible_reactions:
            rxn_hash = CatRxnDescriptorValue.objects.get(descriptor__pk=6793,reaction=rxn).value.value
            score_dict[rxn] = self.hash_MI[rxn_hash] * dist_dict[rxn] #*confidence

        return score_dict


    def create_recommendations(self):
        """Take scored reactions and create RecommendedReaction objects."""
        score_dict = self.score_candidates()

        for rxn, score in score_dict.viewitems():
            print rxn, score
            rec_rxn = RecommendedReaction()
            rec_rxn.reaction_ptr = rxn
            rec_rxn.seed = self.seed
            rec_rxn.nonsense = True # TODO Change in production
            rec_rxn.hidden = True
            rec_rxn.reference = 'Holdingforlaterreplacethis' # What is this supposed to be?
            rec_rxn.labGroup = LabGroup.objects.get(pk=1)
            rec_rxn.score = score
            rec_rxn.saved = True
            rec_rxn.clean()
            rec_rxn.save()

    def recommend(self, k=1000):
        """Recommend k highest scoring plausible reactions."""
        import time
        start_time = time.time()

        self.generator.get_reasonable_compound_amounts(Compound.objects.filter(pk__in=(CompoundQuantity.objects.exclude(amount=None).values_list('compound__pk', flat=True))))
        top_MI_gain_reactant_triples = self.get_best_delta_MI(k=k)

        self.plausible_reactions = self.generator.generate(top_MI_gain_reactant_triples)

        # self.a = top_MI_gain_reactant_triples

        plausible_reaction_ids = [rxn.pk for rxn in self.plausible_reactions][0:100]

        self.plausible_reactions = Reaction.objects.filter(pk__in=plausible_reaction_ids)

        self.sieve = ReactionSieve(self.model_container, self.desired_desc_dict)

        self.plausible_reactions = self.sieve.filter(self.plausible_reactions)

        self.create_recommendations()

        print str(k) + " reactions recommended"
        print("--- %s seconds ---" % (time.time() - start_time))
