from DRP.models import RecommendedReaction
from itertools import product

class ReactionGenerator(object):

    def __init__(self, desc_dict):
      self.desc_dict = desc_dict # A dictionary of (Descriptor, Value List) pairs.

    def generate(self):
      desc_vals = [list(desc_vals.all()) for _, desc_vals in self.desc_dict()]
      # Create a new reaction for each possible pair.
      for new_combo in product(*desc_vals):
        new_rxn = RecommendedReaction()

        # Bind each value to the appropriate descriptor for the new reaction.
        for desc, val in zip(self.desc_dict.keys(), new_combo):
          desc.createValue(new_rxn, val)

        new_rxn.save()
        yield new_rxn

