
class ReactionSieve(object):

  def __init__(self, model_container, desired_desc_dict):
    self.model_container = model_container
    self.desired_desc_dict = desired_desc_dict

  def filter(self, reactions):
    self.model_container.predict(reactions)
    for desc, desc_vals in self.desired_desc_dict.items():
      reactions = filter(lambda rxn: rxn.descriptor_values.filter(value__in=desc_vals,descriptor=desc).exists(), reactions)
    return reactions



