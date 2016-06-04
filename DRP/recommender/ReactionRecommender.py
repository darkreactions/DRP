from DRP.recommender.AbstractRecommender import AbstractRecommender
from DRP.recommender.ReactionGenerator import ReactionGenerator
from DRP.recommender.ReactionSieve import ReactionSieve

# TODO: CompoundQuantities aren't Descriptors, unfortunately.

class Recommender(AbstractRecommender):
  def __init__(self, model_container, desc_dict, desired_desc_dict):
    """
    All Descriptor Dicts use (Descriptor, Value List) pairs, where the
    Value List is a list of any values this Descriptor should have.
    """
    self.generator = ReactionGenerator(desc_dict)
    self.sieve = ReactionSieve(model_container, desired_desc_dict)

  def recommend(self):
    rxns = self.generator.generate()
    desired_rxns = self.sieve.filter(rxns)
    return desired_rxns


"""
# Intended Use:
desc_dict = {DescA:[1,2], DescB:[3], DescC:[4,5,6], OutcomeDesc:[1,2]}
desired_desc_dict = {OutcomeDesc:[2]}
recommender = Recommender(model_container, desc_dict, desired_desc_dict)
good_rxns = recommender.recommend()

"""
