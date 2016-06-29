from copy import deepcopy
from itertools import product

import numpy as np
import xxhash

import DRP
from DRP.plugins.rxndescriptors.rxnhash import calculate
from DRP.models import RecommendedReaction, Compound, CompoundRole, CompoundQuantity, LabGroup, Reaction, NumRxnDescriptorValue


def get_compound_set(role):
    compound_set = set()
    count = 0
    for comp in CompoundQuantity.objects.filter(role=role):
        compound_set.add(comp.compound.id)
        count += 1
        if count > 5:
            break
    return list(compound_set)


class ReactionGenerator(object):
    """Generate the reaction grid space."""

    def __init__(self, desc_dict,  compound_dict=None, amounts_dict=None, compound_roles=None):
        """The desc dict lists descriptors and different values that they may have to construct the grid from"""
        self.desc_dict = desc_dict # A dictionary of (Descriptor, Value List) pairs.
        self.compound_dict = compound_dict
        self.compound_amounts_dict = amounts_dict
        self.compound_roles = compound_roles
        self.triples_dict = {}
        # print self.desc_dict

    def get_all_compound_triples(self):
        org = get_compound_set(6)
        inorg = get_compound_set(7)

        return product(inorg, inorg, org)

    def get_reasonable_compound_amounts(self, compounds):
        """Take a list of compound objects and return a dictionary of compounds to reasonable amounts"""
        compound_dict = {}
        # print list(compounds)
        for compound in compounds:
            quantities = CompoundQuantity.objects.filter(compound=compound).exclude(amount=None)
            quantities = [quantity.amount for quantity in quantities]
            quantities = np.array(quantities)
            mean_quantity = np.mean(quantities)
            quantity_sd = np.std(quantities)
            sample_quantities = [mean_quantity - 2 * quantity_sd, mean_quantity - quantity_sd, mean_quantity,
                                 mean_quantity + quantity_sd, mean_quantity + 2 * quantity_sd]

            sample_quantities = [sq for sq in sample_quantities if sq > 0]
            compound_dict[compound] = sample_quantities

        self.compound_amounts_dict = compound_dict
        # print '~~~~', self.compound_amounts_dict

    def get_all_triple_hashes(self):
        compound_triples = self.get_all_compound_triples()
        triples_dict = {}
        for triple in compound_triples:
            if len(triple) == len(set(triple)):
                h = xxhash.xxh64()  # generates a hash
                compounds = Compound.objects.filter(pk__in=triple).order_by('abbrev')
                for reactant in compounds:
                    h.update(reactant.abbrev)
                triples_dict[h] = compounds
        return triples_dict

    def get_triples_and_amounts(self, compound_triples):
        triples_and_amounts = {}
        # print 'gta', compound_triples
        # print 'cad', self.compound_amounts_dict
        for c1, c2, c3 in compound_triples:
            triples_and_amounts[(c1,c2,c3)] = list(product(self.compound_amounts_dict[c1],self.compound_amounts_dict[c2], self.compound_amounts_dict[c3]))
        
        return triples_and_amounts

    
    def generate(self, compound_triples):
      """Populate the database with every recommended reaction within the grid"""
      # TODO change this line of Casey's code to use .values() method
      #desc_vals = [desc_vals for _, desc_vals in self.desc_dict.iteritems()]
      descs = self.desc_dict.keys()
      desc_vals = self.desc_dict.values()
    
      # Create a new reaction for each possible pair.

      triples_and_amounts = self.get_triples_and_amounts(compound_triples)
      
      desc_combos = product(*desc_vals)
    
      for desc_values_instance in desc_combos:
            for triple in triples_and_amounts:
                for compound_amts in triples_and_amounts[triple]:
                    new_rxn = Reaction()
                    #Currently assigning all of these to Alex's group 
                    new_rxn.labGroup_id = 1 
                    new_rxn.notes = 'Part of a grid search'
                    new_rxn.save()


                    # Generate (numeric) descriptor values for this reaction
                    for desc, dvalue in zip(descs, desc_values_instance):
                        #TODO Make this work for general descriptors
                        new_desc = NumRxnDescriptorValue() 
                        new_desc.id = None
                        new_desc.reaction = new_rxn
                        new_desc.descriptor = desc
                        new_desc.value = dvalue
                        new_desc.save()

                    # Populate the compound quantity fields for the reaction
                    for compound, amount, role in zip(triple, compound_amts, [7,7,6]):
                        compound_quantity = CompoundQuantity()
                        compound_quantity.compound = compound
                        compound_quantity.reaction = new_rxn

                        #TODO Figure out compound role
                        compound_quantity.role = CompoundRole.objects.get(pk=role)
                        #END TODO 
                        compound_quantity.amount = amount

                        compound_quantity.save()

                    new_rxn.save(calcDescriptors=True)

                    # Calculate reaction hash for the proposed reaction
                    calculate(new_rxn)

                    # Save us from calculating the obscene number of reactions in the grid until necessary
                    yield new_rxn
              
        


