from itertools import product

import numpy as np
import xxhash

import DRP
import DRP.plugins.rxndescriptors.drp as drp_rxn_plugin
import DRP.plugins.rxndescriptors.rxnhash as rxnhash_rxn_plugin
from DRP.models import Compound, CompoundRole, CompoundQuantity, Reaction, NumRxnDescriptor, BoolRxnDescriptor, NumRxnDescriptorValue, BoolRxnDescriptorValue, PerformedReaction

INORGANIC_ID = 6
ORGANIC_ID = 7



class ReactionGenerator(object):
    """Generate the reaction grid space."""

    def __init__(self, grid_params, labGroup,  compound_dict=None, amounts_dict=None, compound_roles=None):
        """The desc dict lists descriptors and different values that they may have to construct the grid from"""
        self.grid_params = grid_params # A dictionary of (Descriptor, Value List) pairs.
        self.compound_dict = compound_dict
        self.compound_amounts_dict = amounts_dict
        self.compound_roles = compound_roles
        self.triples_dict = {}
        self.labGroup = labGroup


    def get_compound_set(self, role):
        # TODO: why does this function only select 5?
        compound_set = set()
        count = 0

        valid_performed_reaction_ids = PerformedReaction.objects.filter(valid=True).values_list('id', flat=True).order_by('id')
        reaction_ids = Reaction.objects.filter(labGroup_id=self.labGroup, id__in=valid_performed_reaction_ids).values_list('id', flat=True).order_by('id')
        for comp in CompoundQuantity.objects.filter(role=role, reaction_id__in=reaction_ids):
            print(comp)
            compound_set.add(comp.compound.id)
            if (len(compound_set) > 2):
                break
        print(list(compound_set))
        return list(compound_set)

    def get_all_compound_triples(self):
        # TODO: magic numbers 6 and 7
        org = self.get_compound_set(ORGANIC_ID)
        inorg = self.get_compound_set(INORGANIC_ID)
        # return looks like a tuple of tuple triples
        return product(inorg, inorg, org)

    def get_reasonable_compound_amounts(self, compounds):
        """Take a list of compound objects and return a dictionary of compounds to reasonable amounts"""
        compound_dict = {}
        for compound in compounds:
            quantities = CompoundQuantity.objects.filter(compound=compound).exclude(amount=None)
            quantities = [quantity.amount for quantity in quantities]
            quantities = np.array(quantities)
            mean_quantity = np.mean(quantities)
            quantity_sd = np.std(quantities)
            print("USING ONLY MEAN FOR COMPOUND AMOUNTS CHANGE IN PRODUCTION!")
            sample_quantities = [mean_quantity]
            # TODO: uncomment in production and add more standard deviation slices
            #sample_quantities = [mean_quantity - 2 * quantity_sd, mean_quantity - quantity_sd, mean_quantity,
            #                     mean_quantity + quantity_sd, mean_quantity + 2 * quantity_sd]

            sample_quantities = [sq for sq in sample_quantities if sq > 0]
            compound_dict[compound] = sample_quantities

        # self.compound_amounts_dict is a dictionary with compound objects as keys and a list of values as the value
        self.compound_amounts_dict = compound_dict

    # def get_all_triple_hashes(self):
    #     compound_triples = self.get_all_compound_triples()
    #     triples_dict = {}
    #     for triple in compound_triples:
    #         if len(triple) == len(set(triple)):
    #             h = xxhash.xxh64()  # generates a hash
    #             compounds = Compound.objects.filter(pk__in=triple).order_by('name')
    #             print(compounds)
    #             for reactant in compounds:
    #                 h.update(reactant.name) # problem here!!!, can't resolve keyword 'abbrev' into field, trying just order the reactions by a name
    #             triples_dict[h.hexdigest()] = compounds
    #     print("Triples_dict: ", triples_dict)
    #     # triples_dict has the form of {xxhash: compound triple queryset, ... }
    #     return triples_dict

    def get_all_triples(self):
        compound_triples = self.get_all_compound_triples()
        triples = []
        for triple in compound_triples:
            # make sure our two inorganic reactants aren't the same
            if len(triple) == len(set(triple)):
                # TODO: the ordering here seems to work since triple is already (inorg, inorg, org), but be skeptical
                compounds = Compound.objects.filter(pk__in=triple)
                triples.append(compounds)
        print("Triples: ", triples)
        return triples


    def get_triples_and_amounts(self, compound_triples):
        triples_and_amounts = {}
        for c1, c2, c3 in compound_triples:
            triples_and_amounts[(c1,c2,c3)] = list(product(self.compound_amounts_dict[c1],self.compound_amounts_dict[c2], self.compound_amounts_dict[c3]))
        
        # triples_and_amounts is a dictionary with compound triple tuple as key and tuples with three lists of compound values as the value
        # {(c1,c2,c3): ([value for c1],[value for c2],[value for c3]), ...}
        return triples_and_amounts

    
    def generate(self, compound_triples):
        """Populate the database with every recommended reaction within the grid"""
        grid_params_descs = list(self.grid_params.keys())
        grid_params_desc_vals = list(self.grid_params.values())

        triples_and_amounts = self.get_triples_and_amounts(compound_triples)
        grid_params_desc_combos = product(*grid_params_desc_vals) 
        reaction_set = []
        counter = 0
        """ this big nested loop structure essentially takes each tuple of base descriptor values and each possible triple (and each possible combination of values for said triple)
            and creates a new reaction object, then it loops over the base descriptor values to create the descriptor objects and it loops over the compound values to create compound objects
            finally, after setting up each reaction with it's base descriptor values and compound amounts, we use those to calculate the rest of the descriptors.
        """
        for grid_params_desc_values_instance in grid_params_desc_combos:
            for triple in triples_and_amounts:
                for compound_amts in triples_and_amounts[triple]:
                    new_rxn = Reaction()
                    # TODO: Currently assigning all of these to Alex's group (1) need to change magic number
                    new_rxn.labGroup_id = 1 
                    new_rxn.notes = 'Part of a grid search'
                    new_rxn.save()


                    # Generate (numeric and boolean) base descriptor values for this reaction
                    for desc, dvalue in zip(grid_params_descs, grid_params_desc_values_instance):
                        #TODO Make this work for general descriptors
                        if isinstance(desc, NumRxnDescriptor):
                            new_desc = NumRxnDescriptorValue()
                        if isinstance(desc, BoolRxnDescriptor):
                            new_desc = BoolRxnDescriptorValue()
                        new_desc.id = None
                        new_desc.reaction = new_rxn
                        new_desc.descriptor = desc 
                        new_desc.value = dvalue
                        new_desc.save()

                    # Populate the compound quantity fields for the reaction
                    # Again, magic numbers 6 and 7 that refer to inorganic and organic roles
                    for compound, amount, role in zip(triple, compound_amts, [INORGANIC_ID,INORGANIC_ID,ORGANIC_ID]):
                        compound_quantity = CompoundQuantity()
                        compound_quantity.compound = compound
                        compound_quantity.reaction = new_rxn
                        compound_quantity.role = CompoundRole.objects.get(pk=role)
                        compound_quantity.amount = amount
                        compound_quantity.save()

                    new_rxn.save()
                    counter += 1
                    print('Reaction #{} created.'.format(counter))
                    reaction_set.append(new_rxn)

                    # TODO: Does this need to be a generator since we end up calculating descriptors for all the reactions now...?
                    yield new_rxn
        # generate the rest of the descriptors
        # TODO: find a way to pass a whitelist...
        drp_rxn_plugin.calculate_many(reaction_set)
        rxnhash_rxn_plugin.calculate_many(reaction_set)