import numpy as np
from random import sample

from DRP.models import Compound, CompoundQuantity, CompoundRole, NumericDescriptor

class GridValueGenerator:
    """Currently just sets some expected grid steps as constants and searches for reactants and there ammounts to vary."""
    def __init__(self):
        self.desc_dict = {}
        self.compounddict = {}
        self.compound_amounts_dict = {}
        self.compound_roles = {} #Change if a compound can have more than one role in a reaction
        
    def vary_manual(self):
        self.desc_dict[NumericDescriptor.objects.get(pk=6)] = [1, 2, 3, 8] #, 5, 7, 8]
        self.desc_dict[NumericDescriptor.objects.get(pk=5)] = [1440, 3600]#, 4140, 5040]
        self.desc_dict[NumericDescriptor.objects.get(pk=4)] = [363.15, 423.15]#, 473.15]
        
    def vary_compound(self):
        inorg1 = None
        inorg2 = None
        org = None
        #TODO Add water / solvent
        
        a = Compound.objects.filter(chemicalClasses=4)
        
        org = sample(a,30)
                    
        b = Compound.objects.filter(chemicalClasses=1)
        
        #assuming inorg = inorg2
        inorg = sample(b, 30)

        inorg2= inorg[:15]
        inorg1=inorg[15:]
        
        
        def get_reasonable_compound_amounts(compounds, role):
            """Take a list of compound objects and return a dictionary of compounds to reasonable amounts"""
            compound_dict = {}
            for compound in compounds:
                quantities = CompoundQuantity.objects.filter(compound=compound, role=role).exclude(amount=None)
                self.compound_roles[compound] = CompoundRole.objects.get(pk = role)
                quantities = [quantity.amount for quantity in quantities]
                quantities = np.array(quantities)
#                 try:
                mean_quantity = np.mean(quantities)
#                 except:
#                     for i in quantities:
#                         print i
#                         assert(isinstance(i, float))
#                     assert(False)
                quantity_sd = np.std(quantities)
                sample_quantities = [mean_quantity - 2 * quantity_sd, mean_quantity - quantity_sd, mean_quantity, mean_quantity + quantity_sd, mean_quantity + 2 * quantity_sd]
                sample_quantities = [mean_quantity - quantity_sd]
                sample_quantities = [sq for sq in sample_quantities if sq > 0]
                compound_dict[compound] = sample_quantities
            return compound_dict
        
        self.compounddict["inorg1"] = inorg1
        self.compounddict["inorg2"] = inorg2
        self.compounddict["org"] = org

        
        inorg1 = get_reasonable_compound_amounts(inorg1, 6)
        inorg2 = get_reasonable_compound_amounts(inorg2, 6)
        org = get_reasonable_compound_amounts(org, 7)

        self.compound_amounts_dict.update(inorg1)
        self.compound_amounts_dict.update(inorg2)
        self.compound_amounts_dict.update(org)
        
    def main(self):
        self.vary_manual()
        self.vary_compound()
        return self.desc_dict, self.compounddict, self.compound_amounts_dict, self.compound_roles
        
