from DRP.recommender.ReactionRecommender import Recommender as ReactionRecommender
from DRP.models.rxnDescriptorValues import RxnDescriptorValue

class Recommender(ReactionRecommender):
  def __init__(self, model_container, seed_reaction, desired_desc_dict):
    """
    All Descriptor Dicts use (Descriptor, Value List) pairs, where the
    Value List is a list of any values this Descriptor should have.
    """
    desc_vals = RxnDescriptorValue.objects.filter(reaction=seed_reaction)
    desc_dict = {desc_val.descriptor:desc_val.val for desc_val in desc_vals}
    super(Recommender, self).__init__(model_container, desc_dict, desired_desc_dict)

  def recommend(self):
    return super(Recommender, self).recommend()
