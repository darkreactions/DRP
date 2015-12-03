"""Basic reaction descriptors calculation module"""
import DRP
from itertools import chain
from numpy import mean, average as wmean
from scipy.stats import gmean
from django.db.models import Sum
from utils import setup
import xxhash

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




def calculate(reaction):
    """Calculate the descriptors for this plugin."""

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
                
    
    #Set up the actual descriptor dictionary.
    descriptorDict = setup(_descriptorDict)

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

    for element in elements:
        n = num.objects.get_or_create(
                            reaction=reaction,
                            descriptor=descriptorDict[element + '_mols'],
                            )[0]
        if any(quantity.amount is None for quantity in allCompoundQuantities):
            n.value = None
        else: 
            n.value=sum((quantity.compound.elements[element]['stoichiometry'] * quantity.amount if element in quantity.compound.elements else 0) for quantity in allCompoundQuantities)
        n.save()

    for compoundRole in DRP.models.CompoundRole.objects.all():
        roleQuantities = allCompoundQuantities.filter(role=compoundRole)

        #  number of species in reaction with this role
        n = num.objects.get_or_create(
            reaction=reaction,
            descriptor=descriptorDict['{}_amount_count'.format(compoundRole.label)],
            )[0]
        n.value = roleQuantities.count()
        n.save()

        #  moles ofsum(quantity.amount for quantity in roleQuantities) reactant filling this role in this reaction
        if any(quantity.amount is None for quantity in roleQuantities):
            roleMoles = None
        else:
            roleMoles = sum(quantity.amount for quantity in roleQuantities)
        n= num.objects.get_or_create(
            reaction=reaction,
            descriptor=descriptorDict['{}_amount_molarity'.format(compoundRole.label)],
            )[0]
        n.value=roleMoles
        n.save()

        if roleQuantities.exists():
            for descriptor in DRP.models.NumMolDescriptor.objects.all():
                descriptorValues = DRP.models.NumMolDescriptorValue.objects.filter(compound__in=[quantity.compound for quantity in roleQuantities])
                #  Only do the calculation if the right number of descriptor values are present and all of them are not NULL
                if descriptorValues.count() == roleQuantities.count() and not any(descriptorValue.value is None for descriptorValue in descriptorValues):
                    n = num.objects.get_or_create(
                        reaction=reaction,
                        descriptor=descriptorDict['{}_{}_{}'.format(compoundRole.label, descriptor.csvHeader, 'Max')],
                        
                    )[0]
                    if any(descriptorValue.value is None for descriptorValue in descriptorValues):
                        n.value=None
                    else:
                        n.value=max(descriptorValue.value for descriptorValue in descriptorValues)
                    n.save()
                    num.objects.get_or_create(
                        reaction=reaction,
                        descriptor=descriptorDict['{}_{}_{}'.format(compoundRole.label, descriptor.csvHeader, 'Range')],
                    )[0]
                    if any(descriptorValue.value is None for descriptorValue in descriptorValues):
                        n.value=None
                    else:
                        n.value=max(descriptorValue.value for descriptorValue in descriptorValues) - min(descriptorValue.value for descriptorValue in descriptorValues)
                    n.save()
                    n=num.objects.get_or_create(
                        reaction=reaction,
                        descriptor=descriptorDict['{}_{}_{}_{}'.format(compoundRole.label, descriptor.csvHeader, 'gmean', 'molarity')],
                    )[0]
                    if any(descriptorValues.get(compound=quantity.compound).value is None for quantity in roleQuantities) or roleMoles == 0 or roleMoles is None:
                        n.value=None
                    else:
                        n.value=gmean(list(descriptorValues.get(compound=quantity.compound).value*(quantity.amount/roleMoles) for quantity in roleQuantities))
                    n.save()
                    n=num.objects.get_or_create(
                        reaction=reaction,
                        descriptor=descriptorDict['{}_{}_{}_{}'.format(compoundRole.label, descriptor.csvHeader, 'gmean', 'count')],
                    )[0]
                    n.value=gmean(list(descriptorValues.get(compound=quantity.compound).value for quantity in roleQuantities))
                    n.save()
            for descriptor in DRP.models.OrdMolDescriptor.objects.all():
                descriptorValues = DRP.models.OrdMolDescriptorValue.objects.filter(compound__in=[quantity.compound for quantity in roleQuantities])
                #  Only do the calculation if the right number of descriptor values are present and all of them are not NULL
                if descriptorValues.count() == roleQuantities.count() and not any(descriptorValue.value is None for descriptorValue in descriptorValues):
                    for i in range(descriptor.minimum, descriptor.maximum+1): #  because Still python...
                        n=num.objects.get_or_create(
                            reaction=reaction,
                            descriptor=descriptorDict['{}_{}_{}_count'.format(compoundRole.label, descriptor.csvHeader, i)],
                        )[0]
                        n.value=sum(1 for value in descriptorValues if value.value == i)
                        n.save()
                        n=num.objects.get_or_create(
                            reaction=reaction,
                            descriptor=descriptorDict['{}_{}_{}_molarity'.format(compoundRole.label, descriptor.csvHeader, i)],
                        )[0]
                        quantities = roleQuantities.filter(compound__ordmoldescriptorvalue__value=i, compound__ordmoldescriptorvalue__descriptor__pk=descriptor.pk)
                        if any(quantity.amount is None for quantity in quantities):
                            n.value=None
                        else:
                            n.value=sum(quantity.amount for quantity in quantities)
                        n.save()
            for descriptor in DRP.models.BoolMolDescriptor.objects.all():
                descriptorValues = DRP.models.BoolMolDescriptorValue.objects.filter(compound__in=[quantity.compound for quantity in roleQuantities])
                if descriptorValues.count() == roleQuantities.count() and not any(descriptorValue.value is None for descriptorValue in descriptorValues):
                    for i in (True, False): #  because Still python...
                        n=num.objects.get_or_create(
                            reaction=reaction,
                            descriptor=descriptorDict['{}_{}_{}_count'.format(compoundRole.label, descriptor.csvHeader, i)],
                        )[0]
                        n.value=sum(1 for value in descriptorValues if value.value == i)
                        n.save()
                        n=num.objects.get_or_create(
                            reaction=reaction,
                            descriptor=descriptorDict['{}_{}_{}_molarity'.format(compoundRole.label, descriptor.csvHeader, i)],
                        )[0]
                        quantities = roleQuantities.filter(compound__boolmoldescriptorvalue__value=i, compound__boolmoldescriptorvalue__descriptor__pk=descriptor.pk)
                        if any(quantity.amount is None for quantity in quantities):
                            n.value = None
                        else:
                            n.value=sum(quantity.amount for quantity in quantities)
                        n.save()
            for descriptor in DRP.models.CatMolDescriptor.objects.all():
                descriptorValues = DRP.models.CatMolDescriptorValue.objects.filter(compound__in=[quantity.compound for quantity in roleQuantities])
                if descriptorValues.count() == roleQuantities.count() and not any(descriptorValue.value is None for descriptorValue in descriptorValues):
                    for permValue in descriptor.permittedValues.all():
                        n=num.objects.get_or_create(
                            reaction=reaction,
                            descriptor=descriptorDict['{}_{}_{}_count'.format(compoundRole.label, descriptor.csvHeader, permValue.value)],
                        )[0]
                        n.value=descriptorValues.filter(value=permValue).count()
                        n.save()
                        n=num.objects.get_or_create(
                            reaction=reaction,
                            descriptor=descriptorDict['{}_{}_{}_molarity'.format(compoundRole.label, descriptor.csvHeader, permValue.value)],
                        )[0]
                        quantities = roleQuantities.filter(compound__catmoldescriptorvalue__value=permValue, compound__catmoldescriptorvalue__descriptor__pk=descriptor.pk)
                        if any(quantity.amount is None for quantity in quantities):
                            n.value = None
                        else:
                            n.value=sum(quantity.amount for quantity in quantities)
                        n.save()
