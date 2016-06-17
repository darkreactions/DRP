"""I think this is a first attempt at a crude filter for reactions... PA."""

class ReactionSieve(object):

  """Reaction seive class."""

  def __init__(self, model_container, desired_desc_dict):
    """Initialiser."""
    self.model_container = model_container
    self.desired_desc_dict = desired_desc_dict

  def filter(self, reactions):
    """Filtration function. Uses a lambda function. Suspect."""
    self.model_container.predict(reactions)
    for desc, desc_vals in self.desired_desc_dict.items():
      reactions = filter(lambda rxn: rxn.descriptor_values.filter(value__in=desc_vals, descriptor=desc).exists(), reactions)
    return reactions



