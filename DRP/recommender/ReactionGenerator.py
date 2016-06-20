from copy import deepcopy
from itertools import product

from DRP.models import RecommendedReaction, CompoundQuantity, LabGroup, Reaction, NumRxnDescriptorValue

class ReactionGenerator(object):
    """Generate the reaction grid space."""

    def __init__(self, desc_dict, compound_dict, amounts_dict, compound_roles):
        """The desc dict lists descriptors and different values that they may have to construct the grid from"""
        self.desc_dict = desc_dict # A dictionary of (Descriptor, Value List) pairs.
        self.compound_dict = compound_dict
        self.compound_amounts_dict = amounts_dict
        self.compound_roles = compound_roles


    def get_triples_and_amounts(self):
        lists_of_compounds = self.compound_dict.values()
        compound_triples = product(*lists_of_compounds)
        triples_and_amounts = {}
        for c1, c2, c3 in compound_triples:
            triples_and_amounts[(c1,c2,c3)] =  list(product(self.compound_amounts_dict[c1],self.compound_amounts_dict[c2], self.compound_amounts_dict[c3]))
        
        return triples_and_amounts

    
    def generate(self):
      """Populate the database with every recommended reaction within the grid"""
      # TODO change this line of Casey's code to use .values() method
      #desc_vals = [desc_vals for _, desc_vals in self.desc_dict.iteritems()]
      descs = self.desc_dict.keys()
      desc_vals = self.desc_dict.values()
    
      # Create a new reaction for each possible pair.

      triples_and_amounts = self.get_triples_and_amounts()
      
      desc_combos = product(*desc_vals)
    
      for desc_values_instance in desc_combos:
            for triple in triples_and_amounts:
                for compound_amts in triples_and_amounts[triple]:
                    new_rxn = Reaction()
                    new_rxn.labGroup_id = 1 # LabGroup.objects.get(pk=1)
                    new_rxn.notes = 'Part of a grid search'
                    new_rxn.save()


                    for desc, dvalue in zip(descs, desc_values_instance):
                        #TODO Make this work for general descriptors
                        new_desc = NumRxnDescriptorValue() 
                        new_desc.id = None
                        new_desc.reaction = new_rxn
                        new_desc.descriptor = desc
                        new_desc.value = dvalue
                        new_desc.save()
                   


                    for compound, amount in zip(triple, compound_amts):
                        compound_quantity = CompoundQuantity()
                        compound_quantity.compound = compound
                        compound_quantity.reaction = new_rxn

                        #TODO Figure out compound role
                        compound_quantity.role = self.compound_roles[compound]
                        #END TODO 
                        compound_quantity.amount = amount

                        compound_quantity.save()
                       
                    new_rxn.save(calcDescriptors=True)

                    yield new_rxn
              
        


