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
        xxhash.update(reactant.abbrev)
    p = perm.objects.get_or_create(descriptor=descriptorDict['rxnSpaceHash1'], value=h)[0]
    c = cat.objects.get_or_create(reaction=reaction,descriptor=descriptorDict['rxnSpaceHash1'])[0] 
    c.value = p
    c.save()

    # Calculate the elemental molarities
    allCompoundQuantities = CompoundQuantity.objects.filter(reaction=reaction)

    for element in elements:
        num.objects.get_or_create(
                            reaction=reaction,
                            descriptor=descriptorDict[element + '_mols'],
                            value=sum(quantity.compound.elements[element]['stoichiometry'] * quantity.amount for quantity in allCompoundQuantities))

    for compoundRole in DRP.models.CompoundRole.objects.all():
        roleQuantities = allCompoundQuantities.filter(role=compoundRole)

        #  number of species in reaction with this role
        num.objects.get_or_create(
            reaction=reaction,
            descriptor=descriptorDict['{}_amount_count'.format(compoundRole.label)],
            value=roleQuantities.count())

        #  moles ofsum(quantity.amount for quantity in roleQuantities) reactant filling this role in this reaction
        roleMoles = sum(quantity.amount for quantity in roleQuantities)
        num.objects.get_or_create(
            reaction=reaction,
            descriptor=descriptorDict['{}_amount_molarity'.format(compoundRole.label)],
            value=roleMoles)

        if roleQuantities.exists():
            for descriptor in DRP.models.NumMolDescriptor.objects.all():
                descriptorValues = DRP.models.NumMolDescriptorValue.objects.filter(compound__in=[quantity.compound for quantity in roleQuantities])
                #  Only do the calculation if the right number of descriptor values are present and all of them are not NULL
                if descriptorValues.count() == roleQuantities.count() and not any(descriptorValue.value is None for descriptorValue in descriptorValues):
                    num.objects.get_or_create(
                        reaction=reaction,
                        descriptor=descriptorDict['{}_{}_{}'.format(compoundRole.label, descriptor.csvHeader, 'Max')],
                        value=max(descriptorValue.value for descriptorValue in descriptorValues)
                    )
                    num.objects.get_or_create(
                        reaction=reaction,
                        descriptor=descriptorDict['{}_{}_{}'.format(compoundRole.label, descriptor.csvHeader, 'Range')],
                        value=max(descriptorValue.value for descriptorValue in descriptorValues) - min(descriptorValue.value for descriptorValue in descriptorValues)
                    )
                    num.objects.get_or_create(
                        reaction=reaction,
                        descriptor=descriptorDict['{}_{}_{}_{}'.format(compoundRole.label, descriptor.csvHeader, 'gmean', 'molarity')],
                        value=gmean(descriptorValues.get(compound=quantity.compound).value*(quantity.amount/roleMoles) for quantity in roleQuantities)
                    ) 
                    num.objects.get_or_create(
                        reaction=reaction,
                        descriptor=descriptorDict['{}_{}_{}_{}'.format(compoundRole.label, descriptor.csvHeader, 'gmean', 'count')],
                        value=gmean(descriptorValues.get(compound=quantity.compound).value for quantity in roleQuantities)
                    )
            for descriptor in DRP.models.OrdMolDescriptor.objects.all():
                descriptorValues = DRP.models.OrdMolDescriptorValue.objects.filter(compound__in=[quantity.compound for quantity in roleQuantities])
                #  Only do the calculation if the right number of descriptor values are present and all of them are not NULL
                if descriptorValues.count() == roleQuantities.count() and not any(descriptorValue.value is None for descriptorValue in descriptorValues):
                    for i in range(descriptor.minimum, descriptor.maximum+1): #  because Still python...
                        num.objects.get_or_create(
                            reaction=reaction,
                            descriptor=descriptorDict['{}_{}_{}_count'.format(compoundRole.label, descriptor.csvHeader, i)],
                            value=len(value for value in descriptorValues if value.value == i)
                        )
                        num.objects.get_or_create(
                            reaction=reaction,
                            descriptor=descriptorDict['{}_{}_{}_molarity'.format(compoundRole.label, descriptor.csvHeader, i)],
                            value=sum(quantity.amount for quantity in roleQuantities.filter(compound__ordmoldescriptorvalue__value=i, compound__ordmoldescriptorvalue__descriptor__pk=descriptor.pk))
                        )
            for descriptor in DRP.models.BoolMolDescriptor.objects.all():
                descriptorValues = DRP.models.BoolMolDescriptorValue.objects.filter(compound__in=[quantity.compound for quantity in roleQuantities])
                if descriptorValues.count() == roleQuantities.count() and not any(descriptorValue.value is None for descriptorValue in descriptorValues):
                    for i in (True, False): #  because Still python...
                        num.objects.get_or_create(
                            reaction=reaction,
                            descriptor=descriptorDict['{}_{}_{}_count'.format(compoundRole.label, descriptor.csvHeader, )],
                            value=len(value for value in descriptorValues if value.value == i)
                        )
                        num.objects.get_or_create(
                            reaction=reaction,
                            descriptor=descriptorDict['{}_{}_{}_molarity'.format(compoundRole.label, descriptor.csvHeader, i)],
                            value=sum(quantity.amount for quantity in roleQuantities.filter(compound__boolmoldescriptorvalue__value=i, compound__boolmoldescriptorvalue__descriptor__pk=descriptor.pk))
                        )
            for descriptor in DRP.models.CatMolDescriptor.objects.all():
                descriptorValues = DRP.models.CatMolDescriptorValue.objects.filter(compound__in=[quantity.compound for quantity in roleQuantities])
                if descriptorValues.count() == roleQuantities.count() and not any(descriptorValue.value is None for descriptorValue in descriptorValues):
                    for permValue in descriptor.permittedValues.all():
                        num.objects.get_or_create(
                            reaction=reaction,
                            descriptor=descriptorDict['{}_{}_{}_count'.format(compoundRole.label, descriptor.csvHeader, permValue.value)],
                            value=len(value for value in descriptorValues.filter(value=permValue))
                        )
                        num.objects.get_or_create(
                            reaction=reaction,
                            descriptor=descriptorDict['{}_{}_{}_molarity'.format(compoundRole.label, descriptor.csvHeader, permValue.value)],
                            value=sum(quantity.amount for quantity in roleQuantities.filter(compound__catmoldescriptorvalue__value=permValue, compound__catmoldescriptorvalue__descriptor__pk=descriptor.pk))
                        )
