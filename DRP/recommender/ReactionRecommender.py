from copy import deepcopy

import numpy as np
from sklearn.metrics.cluster import adjusted_mutual_info_score

from DRP.plugins.rxndescriptors.rxnhash import calculate
from DRP.recommender.AbstractRecommender import AbstractRecommender
from DRP.recommender.ReactionGenerator import ReactionGenerator
from DRP.recommender.ReactionSieve import ReactionSieve
from DRP.models import PerformedReaction, BoolRxnDescriptorValue, CatRxnDescriptorValue, CompoundQuantity, Reaction, Compound


class Recommender(AbstractRecommender):
    def __init__(self, model_container, grid_params, desired_desc_dict):
        """
        All Descriptor Dicts use (Descriptor, Value List) pairs, where the
        Value List is a list of any values this Descriptor should have.
        """

        print grid_params
        self.generator = ReactionGenerator(grid_params)

        # self.plausible_reactions]
        self.model_container = model_container
        self.desired_desc_dict = desired_desc_dict

        # react_min = min(reaction_ids)
        # react_max = max(reaction_ids)
        # react_range = (react_min, react_max)
        #
        # self.plausible_reactions =  Reaction.objects.filter(pk__range=react_range)

        # self.plausible_reactions = self.sieve.filter(Reaction.objects.filter(pk__range=react_range))
        # self.committed = []
        
        # TODO remove
        self.desc_dict = grid_params


        # [0]
        # self.compound_dict = grid_params[1]

        # self.react_range = react_range

        self.lshash = None


    def build_hash_outcome_dict(self):
        outcome_dict = {}
        count = 0
        for rxn in PerformedReaction.objects.all(): 
            count +=1
            # if count % 500 == 0:
            #     print count, '/', len(PerformedReaction.objects.all())
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
        hashes, outcomes = zip(*list(outcome_dict.items()))
        
        return adjusted_mutual_info_score(hashes, outcomes)

    def hash_perform_generator(self, compound_dict, anash):
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
        base_dict = self.build_hash_outcome_dict()

        triples_dict =  self.generator.get_all_triple_hashes()

        hash_min_delta_MI = []
        for rxn_hash in triples_dict:
            hash_min_delta_MI.append(self.hash_perform_generator(base_dict, rxn_hash))

        hash_min_delta_MI.sort(key=(lambda x: x[1]), reverse=True)

        greatest_delta_MI_compounds = []
        for a_hash in hash_min_delta_MI[:k]:
            greatest_delta_MI_compounds.append(triples_dict[a_hash[0]])

        return greatest_delta_MI_compounds

    def get_lshash(self):
        from lshash import LSHash
        lsh = None
        headers = []
        for i in self.plausible_reactions.rows(True):
            print i
            print len(i)
            for header in i :
                headers.append(header)
        #
        # lsh = LSHash(24, len(headers))
        # for i in PerformedReaction.objects.all().rows(True):
        #     to__index = []
        #     print type(i)
        #     print i.keys()
        #     for header in headers:
        #         print 'iiiii', [header]
        #         print type(header)
        #         to__index.append(i[header])
        #
        #     print 'to_index', to__index
        #     lsh.index(to__index)
        #
        # self.lshash = lsh
        #
        # return headers

    def get_furthest(self,k=100):
        headers = self.get_lshash()
        for rxn in self.plausible_reactions:
            pass


    def recommend(self):

        self.generator.get_reasonable_compound_amounts(Compound.objects.filter(pk__in=(CompoundQuantity.objects.exclude(amount=None).values_list('compound__pk', flat=True))))
        top_MI_gain_reactant_triples = self.get_best_delta_MI(k=12)

        print top_MI_gain_reactant_triples

        print 'step1'

        self.plausible_reactions = self.generator.generate(top_MI_gain_reactant_triples)

        self.a = top_MI_gain_reactant_triples

        print 'step2'

        count = 0
        plausible_reaction_ids = []
        for i in self.plausible_reactions:
            plausible_reaction_ids.append(i.pk)
            if count > 1000:
                break
            else:
                count += 12

        # plausible_reaction_ids = [rxn.pk for rxn in self.plausible_reactions][0:100]

        self.plausible_reactions = Reaction.objects.filter(pk__in=plausible_reaction_ids)

        print 'step3'

        self.sieve = ReactionSieve(self.model_container, self.desired_desc_dict)

        print 'step4'

        self.plausible_reactions = self.sieve.filter(self.plausible_reactions)

        print 'step5'

        self.get_lshash()

        # self.plausible_reactions = self.get_furthest(k=12)
        #
        # self.commit_reactions_to_database()

        return self.plausible_reactions
        #
        # self.plausible_reactions = self.plausible_reactions.filter(pk__in=top_MI_reaction_pks)
        #
        # self.sieve = ReactionSieve(self.model_container, tuple(), self.desired_desc_dict)
        #
        # self.plausible_reactions = self.sieve.filter(self.plausible_reactions)




