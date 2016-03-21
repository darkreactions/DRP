"""Basic reaction descriptors calculation module"""
import DRP
from itertools import chain
from numpy import mean, average as wmean
from scipy.stats import gmean
from django.db.models import Sum
from utils import setup
import xxhash
import warnings

elements = DRP.chemical_data.elements

_descriptorDict = { 
    'rxnSpaceHash1':
        {
            'type': 'cat',
            'name': 'Hash of reaction reactants to partition reaction space',
            'calculatorSoftware': 'DRP/xxhash',
            'calculatorSoftwareVersion': '0.02/{}'.format(xxhash.VERSION),
            'permittedValues':[]
        }
}

#The following adds descriptors to the dictionary in an automated way to save on voluminous code


for element in elements:
    _descriptorDict[element + '_mols'] = {
            'type':'num',
            'name':'Mols of {} in a reaction.'.format(element),
            'calculatorSoftware':'DRP',
            'calculatorSoftwareVersion': '0.02',
            'maximum': None,
            'minimum': 0
        }

#descriptors for generalised aggregation across compound roles


def make_dict():
    weightings = ('molarity', 'count')
    for compoundRole in DRP.models.CompoundRole.objects.all():
        for w in weightings:
            _descriptorDict['{}_amount_{}'.format(compoundRole.label, w)] = {
                    'type': 'num',
                    'name': 'Amount of compounds in the reaction belonging in the role "{}" weighted by {}'.format(compoundRole.label, w),
                    'calculatorSoftware': 'DRP',
                    'calculatorSoftwareVersion': '0.02',
                    'maximum':None,
                    'minimum': 0
                }
        for descriptor in DRP.models.CatMolDescriptor.objects.all():
            for w in weightings:
                for permValue in descriptor.permittedValues.all():
                    _descriptorDict['{}_{}_{}_{}'.format(compoundRole.label, descriptor.csvHeader, permValue.value, w)] = {
                        'type': 'num',
                        'name': 'Amount of reactants in category {} for descriptor "{}" weighted by reactant {} in compound role {}.'.format(
                                permValue.value, descriptor.name, w, compoundRole.label),
                        'calculatorSoftware': 'DRP',
                        'calculatorSoftwareVersion': '0.02',
                        'maximum': None,
                        'minimum': None
                        }
        for descriptor in DRP.models.OrdMolDescriptor.objects.all():
            for w in weightings:
                for i in range(descriptor.minimum, descriptor.maximum+1): #  because python...
                    _descriptorDict['{}_{}_{}_{}'.format(compoundRole.label, descriptor.csvHeader, i, w)] = {
                            'type': 'num',
                            'name': 'Amount of reactants with value {} for descriptor "{}" weighted by reactant {} in compound role {}.'.format(
                                    i, descriptor.name, w, compoundRole.label),
                            'calculatorSoftware': 'DRP',
                            'calculatorSoftwareVersion': '0.02',
                            'maximum': None,
                            'minimum': None
                        }
            for descriptor in DRP.models.BoolMolDescriptor.objects.all():
                for w in weightings:
                    for value in ('True', 'False'):
                        _descriptorDict['{}_{}_{}_{}'.format(compoundRole.label, descriptor.csvHeader, value, w)] = {
                                'type': 'num',
                                'name': 'Amount of reactants with value {} for descriptor "{}" weighted by reactant {} in compound role {}.'.format(
                                        value, descriptor.name, w, compoundRole.label),
                                'calculatorSoftware': 'DRP',
                                'calculatorSoftwareVersion': '0.02',
                                'maximum': None,
                                'minimum': None,
                            }
            for descriptor in DRP.models.NumMolDescriptor.objects.all():
                _descriptorDict['{}_{}_{}'.format(compoundRole.label, descriptor.csvHeader, 'Max')] = {
                    'type': 'num',
                    'name': 'Maximum value for {} aggregated across compounds in role "{}"'.format(descriptor.name, compoundRole.label),
                    'calculatorSoftware': 'DRP',
                    'calculatorSoftwareVersion': 'DRP',
                    'maximum': None,
                    'minimum': None
                    }
                _descriptorDict['{}_{}_{}'.format(compoundRole.label, descriptor.csvHeader, 'Range')] = {
                    'type': 'num',
                    'name': 'Range for {} aggregated across compounds in role "{}"'.format(descriptor.name, compoundRole.label),
                    'calculatorSoftware': 'DRP',
                    'calculatorSoftwareVersion': 'DRP',
                    'maximum': None,
                    'minimum': None
                    }
                for w in weightings:
                    _descriptorDict['{}_{}_{}_{}'.format(compoundRole.label, descriptor.csvHeader, 'gmean', w)] = {
                        'type': 'num',
                        'name': 'Geometric Mean for {} aggregated across compounds in role "{}" normalised by {}'.format(descriptor.name, compoundRole.label, w),
                        'calculatorSoftware': 'DRP',
                        'calculatorSoftwareVersion': 'DRP',
                        'maximum': None,
                        'minimum': None
                        }
    descriptorDict = setup(_descriptorDict)
    return descriptorDict


def delete_descriptors_many(reaction_set, descriptorDict):
    # This could be the same as delete descriptors if we're ok with deleting the descriptor even if not
    # descriptorValues.count() == roleQuantities.count() and not any(descriptorValue.value is None for descriptorValue in descriptorValues)
    # I think yes, but wasn't completely sure, so made the separate function

    descriptors_to_delete = []
    for element in elements:
        descriptors_to_delete.append(descriptorDict[element + '_mols'])


    allCompoundQuantities = DRP.models.CompoundQuantity.objects.filter(reaction__in=reaction_set)
    for compoundRole in DRP.models.CompoundRole.objects.all():
        roleQuantities = allCompoundQuantities.filter(role=compoundRole)
        descriptors_to_delete.append(descriptorDict['{}_amount_count'.format(compoundRole.label)])
        #  number of species in reaction with this role
        descriptors_to_delete.append(descriptorDict['{}_amount_molarity'.format(compoundRole.label)])
        
        if roleQuantities.exists():
            descriptorValues = DRP.models.NumMolDescriptorValue.objects.filter(compound__in=[quantity.compound for quantity in roleQuantities])
            for descriptor in DRP.models.NumMolDescriptor.objects.all():
                descriptors_to_delete.append(descriptorDict['{}_{}_{}'.format(compoundRole.label, descriptor.csvHeader, 'Max')])
                descriptors_to_delete.append(descriptorDict['{}_{}_{}'.format(compoundRole.label, descriptor.csvHeader, 'Range')])
                descriptors_to_delete.append(descriptorDict['{}_{}_{}_{}'.format(compoundRole.label, descriptor.csvHeader, 'gmean', 'molarity')])
                descriptors_to_delete.append(descriptorDict['{}_{}_{}_{}'.format(compoundRole.label, descriptor.csvHeader, 'gmean', 'count')])
            for descriptor in DRP.models.OrdMolDescriptor.objects.all():
                for i in range(descriptor.minimum, descriptor.maximum+1): #  because Still python...
                    descriptors_to_delete.append(descriptorDict['{}_{}_{}_count'.format(compoundRole.label, descriptor.csvHeader, i)])
                    descriptors_to_delete.append(descriptorDict['{}_{}_{}_molarity'.format(compoundRole.label, descriptor.csvHeader, i)])
            for descriptor in DRP.models.BoolMolDescriptor.objects.all():
                for i in (True, False): #  because Still python...
                    descriptors_to_delete.append(descriptorDict['{}_{}_{}_count'.format(compoundRole.label, descriptor.csvHeader, i)])
                    descriptors_to_delete.append(descriptorDict['{}_{}_{}_molarity'.format(compoundRole.label, descriptor.csvHeader, i)])
            for descriptor in DRP.models.CatMolDescriptor.objects.all():
                for permValue in descriptor.permittedValues.all():
                    descriptors_to_delete.append(descriptorDict['{}_{}_{}_count'.format(compoundRole.label, descriptor.csvHeader, permValue.value)])
                    descriptors_to_delete.append(descriptorDict['{}_{}_{}_molarity'.format(compoundRole.label, descriptor.csvHeader, permValue.value)])

    DRP.models.NumRxnDescriptorValue.objects.filter(reaction__in=reaction_set, descriptor__in=descriptors_to_delete).delete()

def delete_descriptors(reaction, descriptorDict):
    allCompoundQuantities = DRP.models.CompoundQuantity.objects.filter(reaction=reaction)

    descriptors_to_delete = []
    for element in elements:
        descriptors_to_delete.append(descriptorDict[element + '_mols'])

    for compoundRole in DRP.models.CompoundRole.objects.all():
        roleQuantities = allCompoundQuantities.filter(role=compoundRole)
        descriptors_to_delete.append(descriptorDict['{}_amount_count'.format(compoundRole.label)])
        #  number of species in reaction with this role
        descriptors_to_delete.append(descriptorDict['{}_amount_molarity'.format(compoundRole.label)])
        
        if roleQuantities.exists():
            for descriptor in DRP.models.NumMolDescriptor.objects.all():
                descriptorValues = DRP.models.NumMolDescriptorValue.objects.filter(compound__in=[quantity.compound for quantity in roleQuantities], descriptor=descriptor)
                #  Only do the calculation if the right number of descriptor values are present and all of them are not NULL
                if descriptorValues.count() == roleQuantities.count() and not any(descriptorValue.value is None for descriptorValue in descriptorValues):
                    descriptors_to_delete.append(descriptorDict['{}_{}_{}'.format(compoundRole.label, descriptor.csvHeader, 'Max')])
                    descriptors_to_delete.append(descriptorDict['{}_{}_{}'.format(compoundRole.label, descriptor.csvHeader, 'Range')])
                    descriptors_to_delete.append(descriptorDict['{}_{}_{}_{}'.format(compoundRole.label, descriptor.csvHeader, 'gmean', 'molarity')])
                    descriptors_to_delete.append(descriptorDict['{}_{}_{}_{}'.format(compoundRole.label, descriptor.csvHeader, 'gmean', 'count')])
            for descriptor in DRP.models.OrdMolDescriptor.objects.all():
                descriptorValues = DRP.models.OrdMolDescriptorValue.objects.filter(compound__in=[quantity.compound for quantity in roleQuantities], descriptor=descriptor)
                #  Only do the calculation if the right number of descriptor values are present and all of them are not NULL
                if descriptorValues.count() == roleQuantities.count() and not any(descriptorValue.value is None for descriptorValue in descriptorValues):
                    for i in range(descriptor.minimum, descriptor.maximum+1): #  because Still python...
                        descriptors_to_delete.append(descriptorDict['{}_{}_{}_count'.format(compoundRole.label, descriptor.csvHeader, i)])
                        descriptors_to_delete.append(descriptorDict['{}_{}_{}_molarity'.format(compoundRole.label, descriptor.csvHeader, i)])
            for descriptor in DRP.models.BoolMolDescriptor.objects.all():
                descriptorValues = DRP.models.BoolMolDescriptorValue.objects.filter(compound__in=[quantity.compound for quantity in roleQuantities], descriptor=descriptor)
                if descriptorValues.count() == roleQuantities.count() and not any(descriptorValue.value is None for descriptorValue in descriptorValues):
                    for i in (True, False): #  because Still python...
                        descriptors_to_delete.append(descriptorDict['{}_{}_{}_count'.format(compoundRole.label, descriptor.csvHeader, i)])
                        descriptors_to_delete.append(descriptorDict['{}_{}_{}_molarity'.format(compoundRole.label, descriptor.csvHeader, i)])
            for descriptor in DRP.models.CatMolDescriptor.objects.all():
                descriptorValues = DRP.models.CatMolDescriptorValue.objects.filter(compound__in=[quantity.compound for quantity in roleQuantities], descriptor=descriptor)
                if descriptorValues.count() == roleQuantities.count() and not any(descriptorValue.value is None for descriptorValue in descriptorValues):
                    for permValue in descriptor.permittedValues.all():
                        descriptors_to_delete.append(descriptorDict['{}_{}_{}_count'.format(compoundRole.label, descriptor.csvHeader, permValue.value)])
                        descriptors_to_delete.append(descriptorDict['{}_{}_{}_molarity'.format(compoundRole.label, descriptor.csvHeader, permValue.value)])

            
    DRP.models.NumRxnDescriptorValue.objects.filter(reaction=reaction, descriptor__in=descriptors_to_delete).delete()



def calculate_many(reaction_set, verbose=False):
    descriptorDict = make_dict()
    if verbose:
        print "Deleting old descriptor values"
    delete_descriptors_many(reaction_set, descriptorDict)
    for i, reaction in enumerate(reaction_set):
        if verbose:
            print "Calculating {} ({}/{})".format(reaction, i+1, len(reaction_set))
        _calculate(reaction, descriptorDict, verbose=verbose)


def calculate(reaction):
    """Calculate the descriptors for this plugin."""
    #Set up the actual descriptor dictionary.
    descriptorDict = make_dict()
    delete_descriptors(reaction, descriptorDict)
    _calculate(reaction, descriptorDict)

def _calculate(reaction, descriptorDict, verbose=False):
    """
    Calculates with the descriptorDict already created
    """

    #descriptor Value classes
    CompoundQuantity = DRP.models.CompoundQuantity
    num = DRP.models.NumRxnDescriptorValue
    cat = DRP.models.CatRxnDescriptorValue
    perm = DRP.models.CategoricalDescriptorPermittedValue
    
    #reaction space descriptor
    h = xxhash.xxh64() #generates a hash
    for reactant in reaction.compounds.all():
        h.update(reactant.abbrev)
    p = perm.objects.get_or_create(descriptor=descriptorDict['rxnSpaceHash1'], value=h.hexdigest())[0]
    c = cat.objects.get_or_create(reaction=reaction,descriptor=descriptorDict['rxnSpaceHash1'])[0] 
    c.value = p
    c.save()

    # Calculate the elemental molarities
    allCompoundQuantities = CompoundQuantity.objects.filter(reaction=reaction)

    vals_to_create = []
    for element in elements:
        n = num(
            reaction=reaction,
            descriptor=descriptorDict[element + '_mols'],
            )
        if any(quantity.amount is None for quantity in allCompoundQuantities):
            n.value = None
        else: 
            n.value=float(sum((quantity.compound.elements[element]['stoichiometry'] * quantity.amount if element in quantity.compound.elements else 0) for quantity in allCompoundQuantities))

        vals_to_create.append(n)


    for compoundRole in DRP.models.CompoundRole.objects.all():
        roleQuantities = allCompoundQuantities.filter(role=compoundRole)
        #  number of species in reaction with this role
        n = num(
            reaction=reaction,
            descriptor=descriptorDict['{}_amount_count'.format(compoundRole.label)],
            )
        n.value = roleQuantities.count()
        vals_to_create.append(n)

        #  moles of sum(quantity.amount for quantity in roleQuantities) reactant filling this role in this reaction
        if any(quantity.amount is None for quantity in roleQuantities):
            roleMoles = None
        else:
            roleMoles = sum(quantity.amount for quantity in roleQuantities)
        n = num(
            reaction=reaction,
            descriptor=descriptorDict['{}_amount_molarity'.format(compoundRole.label)],
            )
        n.value=roleMoles
        vals_to_create.append(n)

        if roleQuantities.exists():
            for descriptor in DRP.models.NumMolDescriptor.objects.all():
                descriptorValues = DRP.models.NumMolDescriptorValue.objects.filter(compound__in=[quantity.compound for quantity in roleQuantities], descriptor=descriptor)
                #  Only do the calculation if the right number of descriptor values are present and all of them are not NULL
                # TODO XXX this count now takes longer than actually creating the value. Can we remove it? (Same below)
                # Ah looks like this is silently skipping inorganic properties for organics. That seems like a bad way to do it...
                # Means I have to silence the warnings below (which previously lead me to a whole bunch of uncalculated stuff.
                # I really don't like failing silently and it seems to have gotten us into some trouble
                if descriptorValues.count() == roleQuantities.count() and not any(descriptorValue.value is None for descriptorValue in descriptorValues):
                    n = num(
                        reaction=reaction,
                        descriptor=descriptorDict['{}_{}_{}'.format(compoundRole.label, descriptor.csvHeader, 'Max')],
                    )
                    if any(descriptorValue.value is None for descriptorValue in descriptorValues):
                        n.value = None
                    else:
                        n.value = max(descriptorValue.value for descriptorValue in descriptorValues)
                    vals_to_create.append(n)
                    n = num(
                        reaction=reaction,
                        descriptor=descriptorDict['{}_{}_{}'.format(compoundRole.label, descriptor.csvHeader, 'Range')],
                    )
                    if any(descriptorValue.value is None for descriptorValue in descriptorValues):
                        n.value=None
                    else:
                        n.value=max(descriptorValue.value for descriptorValue in descriptorValues) - min(descriptorValue.value for descriptorValue in descriptorValues)
                    vals_to_create.append(n)
                    n = num(
                        reaction=reaction,
                        descriptor=descriptorDict['{}_{}_{}_{}'.format(compoundRole.label, descriptor.csvHeader, 'gmean', 'molarity')],
                    )
                    if any(descriptorValues.get(compound = quantity.compound).value is None for quantity in roleQuantities) or roleMoles == 0 or roleMoles is None:
                        n.value=None
                    else:
                        n.value = gmean(list(descriptorValues.get(compound=quantity.compound).value*(quantity.amount/roleMoles) for quantity in roleQuantities))
                    vals_to_create.append(n)
                    n = num(
                        reaction=reaction,
                        descriptor=descriptorDict['{}_{}_{}_{}'.format(compoundRole.label, descriptor.csvHeader, 'gmean', 'count')],
                    )
                    n.value = gmean(list(descriptorValues.get(compound=quantity.compound).value for quantity in roleQuantities))
                    vals_to_create.append(n)
                #elif descriptorValues.count() != roleQuantities.count():
                    #warnings.warn("Skipping {} because there are {} descriptorValues and {} roleQuantities".format(descriptor.heading, descriptorValues.count(), roleQuantities.count()))
                #else:
                    #warnings.warn("Skipping {} because some descriptorValues are None".format(descriptor.heading))
            for descriptor in DRP.models.OrdMolDescriptor.objects.all():
                #  Only do the calculation if the right number of descriptor values are present and all of them are not NULL
                # TODO XXX this count now takes longer than actually creating the value. Can we remove it? (Same below)
                descriptorValues = DRP.models.OrdMolDescriptorValue.objects.filter(compound__in=[quantity.compound for quantity in roleQuantities], descriptor=descriptor)
                if descriptorValues.count() == roleQuantities.count() and not any(descriptorValue.value is None for descriptorValue in descriptorValues):
                    for i in range(descriptor.minimum, descriptor.maximum+1): #  because Still python...
                        n = num(
                            reaction=reaction,
                            descriptor=descriptorDict['{}_{}_{}_count'.format(compoundRole.label, descriptor.csvHeader, i)],
                        )
                        n.value = sum(1 for value in descriptorValues if value.value == i)
                        vals_to_create.append(n)
                        n = num(
                            reaction=reaction,
                            descriptor=descriptorDict['{}_{}_{}_molarity'.format(compoundRole.label, descriptor.csvHeader, i)],
                        )
                        quantities = roleQuantities.filter(compound__ordmoldescriptorvalue__value=i, compound__ordmoldescriptorvalue__descriptor__pk=descriptor.pk)
                        if any(quantity.amount is None for quantity in quantities):
                            n.value=None
                        else:
                            n.value=sum(quantity.amount for quantity in quantities)
                        vals_to_create.append(n)
            for descriptor in DRP.models.BoolMolDescriptor.objects.all():
                descriptorValues = DRP.models.BoolMolDescriptorValue.objects.filter(compound__in=[quantity.compound for quantity in roleQuantities], descriptor=descriptor)
                if descriptorValues.count() == roleQuantities.count() and not any(descriptorValue.value is None for descriptorValue in descriptorValues):
                    for i in (True, False): #  because Still python...
                        n = num(
                            reaction=reaction,
                            descriptor=descriptorDict['{}_{}_{}_count'.format(compoundRole.label, descriptor.csvHeader, i)],
                        )
                        n.value=sum(1 for value in descriptorValues if value.value == i)
                        vals_to_create.append(n)
                        n = num(
                            reaction=reaction,
                            descriptor=descriptorDict['{}_{}_{}_molarity'.format(compoundRole.label, descriptor.csvHeader, i)],
                        )
                        quantities = roleQuantities.filter(compound__boolmoldescriptorvalue__value=i, compound__boolmoldescriptorvalue__descriptor__pk=descriptor.pk)
                        if any(quantity.amount is None for quantity in quantities):
                            n.value = None
                        else:
                            n.value=sum(quantity.amount for quantity in quantities)
                        vals_to_create.append(n)
            for descriptor in DRP.models.CatMolDescriptor.objects.all():
                descriptorValues = DRP.models.CatMolDescriptorValue.objects.filter(compound__in=[quantity.compound for quantity in roleQuantities], descriptor=descriptor)
                if descriptorValues.count() == roleQuantities.count() and not any(descriptorValue.value is None for descriptorValue in descriptorValues):
                    for permValue in descriptor.permittedValues.all():
                        n = num(
                            reaction=reaction,
                            descriptor=descriptorDict['{}_{}_{}_count'.format(compoundRole.label, descriptor.csvHeader, permValue.value)],
                        )
                        n.value=descriptorValues.filter(value=permValue).count()
                        vals_to_create.append(n)
                        n = num(
                            reaction=reaction,
                            descriptor=descriptorDict['{}_{}_{}_molarity'.format(compoundRole.label, descriptor.csvHeader, permValue.value)],
                        )
                        quantities = roleQuantities.filter(compound__catmoldescriptorvalue__value=permValue, compound__catmoldescriptorvalue__descriptor__pk=descriptor.pk)
                        if any(quantity.amount is None for quantity in quantities):
                            n.value = None
                        else:
                            n.value = sum(quantity.amount for quantity in quantities)
                        vals_to_create.append(n)
    if verbose:
        print "Creating {} values".format(len(vals_to_create))
    num.objects.bulk_create(vals_to_create)
