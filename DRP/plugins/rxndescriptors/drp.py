"""Basic reaction descriptors calculation module"""
from itertools import chain
from numpy import mean, average as wmean
from scipy.stats import gmean
from django.conf import settings
from django.models import Sum
from utils import setup
import xxhash
import DRP 

elements = DRP.chemical_data.elements

_descriptorDict = { 
    'numberInorganic':
        {
            'type': 'num',
            'name': 'Number of Inorganic Components',
            'calculatorSoftware': 'DRP',
            'calculatorSoftwareVersion': '0.02',
            'maximum': None,
            'minimum': 0,
        },
    'inorgWaterMolRatio':
        {
            'type':'num',
            'name':'Inorganic:Water molar ratio',
            'calculatorSoftware':'DRP',
            'calculatorSoftwareVersion': '0.02',
            'maximum': None,
            'minimum': 0
        },
    'inorgOrgMolRatio':
        {
            'type': 'num',
            'name': 'Inorganic:Organic molar ratio',
            'calculatorSoftware': 'DRP',
            'calculatorSoftwareVersion': '0.02',
            'maximum': None,
            'minimum': 0
        },
    'notWaterWaterMolRatio':
        {
            'type': 'num',
            'name': 'Water:Not Water molar ratio',
            'calculatorSoftware': 'DRP',
            'calculatorSoftwareVersion': '0.02',
            'maximum': None,
            'minimum': 0
        }
    'rxnSpaceHash1':
        {
            'type': 'cat',
            'name': 'Hash of reaction reactants to partition reaction space',
            'calculatorSoftware': 'DRP/xxhash',
            'calculatorSoftwareVersion': '0.02/{}'.format(xxhash.VERSION)
        }
}

#The following adds descriptors to the dictionary in an automated way to save on voluminous code

#Inorganic descriptors
inorgAtomicProperties = (
    'ionization_energy',
    'electron_affinity',
    'pauling_electronegativity',
    'pearson_electronegativity',
    'hardness',
    'atomic_radius'
)
weightings = (
    ('', ''), 
    ('MolWeighted', ' weighted by the molarity of the relevant atoms'),
    ('StoichWeighted', ' weighted in each species stoichiometrically, and aggregated uniformly')
)
aggregations = (
    ('Max', 'Maximum'), 
    ('Range', 'Range'),
    ('Geom', 'Geometric Average')
)


for prop in inorgAtomicProperties:
    stem = 'drpInorgAtom' + prop.title().replace('_', '') 
    for weighting in weightings:
        for agg in aggregations:
             _descriptorDict[stem + weighting[0] + agg[0]] = {
                'type':'num',
                'name': '{} {}{}.'.format(agg[1], prop.replace('_', ' '), weighting[1]),
                'calculatorSoftware':'DRP',
                'calculatorSoftwareVersion':'0.02',
                'maximum':None,
                'minimum':None
            }

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

normalisations = (('molarity', 'mol_norm'), ('count','base_norm'))
for compoundRole in DRP.models.CompoundRole.objects.all():
    for descriptor in DRP.models.CatMolDescriptor.objects.all():
        for norm in normalisations:
            _descriptorDict['{}_{}_{}_{}'.format(compoundRole.label, descriptor.csvHeader, 'mode', norm[1])] = {
                    'type': 'cat',
                    'name': descriptor.name + ' Modal value for "{}" role for reaction normalised by {}.'.format(compoundRole.label, norm[0]),
                    'calculatorSoftware': 'DRP',
                    'calculatorSoftwareVersion': '0.02',
                    'permittedValues': [value.value for value in descriptor.permittedValues.all()]
                }
            for permValue in descriptor.permittedValues.all():
                _descriptorDict['{}_{}_{}_{}'.format(compoundRole.label, descriptor.csvHeader, permValue.value, norm[1])] = {
                        'type': 'num',
                        'name': 'Proportion of reactants in category {} for descriptor "{}" weighted by normalised reactant {} in compound role {}.'.format(
                                permValue.value, descriptor.name, norm[0], compoundRole.label),
                        'calculatorSoftware': 'DRP',
                        'calculatorSoftwareVersion': '0.02',
                        'maximum': 1,
                        'minimum': 0
                    }
    for descriptor in DRP.models.OrdMolDescriptor.objects.all():
        for norm in normalisations:
            _descriptorDict['{}_{}_{}_{}'.format(compoundRole.label, descriptor.csvHeader, 'mode', norm[1]] = {
                    'type': 'ord',
                    'name': descriptor.name + ' Modal value for "{}" role for reaction normalised by {}.'.format(compoundRole.label, norm[0]),
                    'calculatorSoftware': 'DRP',
                    'calculatorSoftwareVersion': '0.02',
                    'maximum': descriptor.maximum,
                    'minimum': descriptor.minimum
                }
            _descriptorDict['{}_{}_{}_{}'.format(compoundRole.label, descriptor.csvHeader, '75pc', norm[1]] {
                    'type': 'ord',
                    'name': descriptor.name + ' 75th percentile value for "{}" role for reaction normalised by {}.'.format(compoundRole.label, norm[0]),
                    'calculatorSoftware': 'DRP',
                    'calculatorSoftwareVersion': '0.02',
                    'maximum': descriptor.maximum,
                    'minimum': descriptor.minimum
                }
            _descriptorDict['{}_{}_{}_{}'.format(compoundRole.label, descriptor.csvHeader, '50pc', norm[1]] {
                    'type': 'ord',
                    'name': descriptor.name + ' 50th percentile value for "{}" role for reaction normalised by {}.'.format(compoundRole.label, norm[0]),
                    'calculatorSoftware': 'DRP',
                    'calculatorSoftwareVersion': '0.02',
                    'maximum': descriptor.maximum,
                    'minimum': descriptor.minimum
                }
            _descriptorDict['{}_{}_{}_{}'.format(compoundRole.label, descriptor.csvHeader, '7525_IPR', norm[1]] {
                    'type': 'ord',
                    'name': descriptor.name + ' interpercentile range for 25th and 75th percentile for "{}" role for reaction normalised by {}.'.format(compoundRole.label, norm[0]),
                    'calculatorSoftware': 'DRP',
                    'calculatorSoftwareVersion': '0.02',
                    'maximum': descriptor.maximum-descriptor.minimum,
                    'minimum': 0
                }
            for i in range(descriptor.minimum, descriptor.maximum+1): #  because python...
                _descriptorDict['{}_{}_{}_{}'.format(compoundRole.label, descriptor.csvHeader, i, norm[1])] = {
                        'type': 'num',
                        'name': 'Proportion of reactants with value {} for descriptor "{}" weighted by normalised reactant {} in compound role {}.'.format(
                                i, descriptor.name, norm[0], compoundRole.label),
                        'calculatorSoftware': 'DRP',
                        'calculatorSoftwareVersion': '0.02',
                        'maximum': 1,
                        'minimum': 0
                    }

        for descriptor in DRP.models.BoolMolDescriptor.objects.all():
            for norm in normalisations:
                for value in ('True', 'False'):
                    _descriptorDict['{}_{}_{}_{}'.format(compoundRole.label, descriptor.csvHeader, value, norm[1])] = {
                            'type': 'num',
                            'name': 'Proportion of reactants with value {} for descriptor "{}" weighted by normalised reactant {} in compound role {}.'.format(
                                    value, descriptor.name, norm[0], compoundRole.label),
                            'calculatorSoftware': 'DRP',
                            'calculatorSoftwareVersion': '0.02',
                            'maximum': 1,
                            'minimum': 0
                        }
                _descriptorDict['{}_{}_{}_{}'.format(compoundRole.label, descriptor.csvHeader, 'mode', norm[1]] = {
                        'type': 'ord',
                        'name': descriptor.name + ' Modal value for "{}" role for reaction normalised by {}.'.format(compoundRole.label, norm[0]),
                        'calculatorSoftware': 'DRP',
                        'calculatorSoftwareVersion': '0.02',
                    }
            _descriptorDict['{}_{}_{}'.format(compoundRole.label, descriptor.csvHeader, 'both_present')] = {
                    'type': 'bool',
                    'name': 'Both true and false presence test for "{}" role for reaction.'.format(compoundRole.label, norm[0]),
                    'calculatorSoftware': 'DRP',
                    'calculatorSoftwareVersion': '0.02'
                }
            

#Set up the actual descriptor dictionary.
descriptorDict = setup(_descriptorDict)

for element, info in elements.items():
    if (element == 'Se') or (info['group'] in range(3, 13)) or ((info['group'] > 12) and ((not info['nonmetal']) or info['metalloid'])):
        inorgElements[element] = info 


def calculate(reaction):
    """Calculate the descriptors for this plugin."""

    #descriptor Value classes
    num = DRP.models.NumRxnDescriptorValue
    cat = DRP.models.CatRxnDescriptorValue
    perm = DRP.models.CategoricalDescriptorPermittedValue

    #reaction space descriptor
    h = xxhash.xxh64() #generates a hash
    for reactant in reaction.compounds:
        xxhash.update(reactant.abbrev)
    p = perm.objects.get_or_create(descriptor=descriptorDict['rxnSpaceHash1'], value=h)
    cat.objects.get_or_create(reaction=reaction,descriptor=descriptorDict['rxnSpaceHash1'], value=p) 

    #inorganic descriptor calculations
    inorgCompoundQuantities = DRP.models.CompoundQuantity.objects.filter(reaction=reaction, role__label='Inorg')  # These two lines get the inorganic compounds for htis reaction.
    sumInorgAmount = inorgCompoundQuantities.aggregate(Sum('amount'))


    MAX = 1
    RANGE = 2
    GMEAN = 3
    WMOL = 1
    WSTOICH = 2

    def inorgAtomicAggregate(propertyLabel, function, weighting=None):
        """Figure out the aggregation for inorganic atom-based properties. This is done a lot so this saves on code repetition."""
        if weighting is None:
            if function == MAX:
                return max(max(inorgElements[element][propertyLabel] for element in quantity.compound.elements if element in inorgElements) for quantity in inorgCompoundQuantities)
            elif function == RANGE:
                maximum = max(max(inorgElements[element][propertyLabel] for element in quantity.compound.elements if element in inorgElements) for quantity in inorgCompoundQuantities)
                minimum = min(min(inorgElements[element][propertyLabel] for element in quantity.compound.elements) if element in inorgElements for quantity in inorgCompoundQuantities)
                return maximum - minimum
            elif function == GMEAN:
                return gmean(chain(inorgElements[element][propertyLabel] for element in quantity.compound.elements if element in inorgElements) for quantity in inorgCompoundQuantities)
            else:
                raise TypeError('Unrecognised function selection constant')
        elif weighting == WMOL:
            if function == MAX:
                return max(max(inorgElements[element][propertyLabel]*(quantity.amount/sumInorgAmount)*info['stoichiometry'] for element, info in quantity.compound.elements if element in inorgElements) for quantity in inorgCompoundQuantities)
            elif function == RANGE:
                maximum = max(max(inorgElements[element][propertyLabel]*(quantity.amount/sumInorgAmount)*info['stoichiometry'] for element, info in quantity.compound.elements if element in inorgElements) for quantity in inorgCompoundQuantities)
                minimum = min(min(inorgElements[element][propertyLabel]*(quantity.amount/sumInorgAmount)*info['stoichiometry'] for element, info in quantity.compound.elements if element in inorgElements) for quantity in inorgCompoundQuantities)
                return maximum-minimum
             elif function == GMEAN:
                return gmean(chain(inorgElements[element][propertyLabel]*(quantity.amount/sumInorgAmount)*info['stoichiometry'] for element, info in quantity.compound.elements if element in inorgElements) for quantity in inorgCompoundQuantities))
            else:
                raise TypeError('Unrecognised function selection constant')
        elif weighting == WSTOICH: 
            values = (wmean((inorgElements[element][propertyLabel] for element in quantity.compound.elements if element in inorgElements), 
                                weights=(info['stoichiometry'] for element, info in quantity.compound.elements.items()))
                                            for quantity in inorgCompoundQuantities)
            if function == MAX:
                return max(values)
            elif function == RANGE:
                return max(values) - min(values)
            elif function == GMEAN:
                return gmean(values)
            else:
                raise TypeError('Unrecognised function selection constant')
        else:
            raise TypeError('Unrecognised weighting scheme')

    # Number of inorganic components
    num.objects.get_or_create(
                            reaction=reaction,
                            descriptor=descriptorDict['numberInorganic'],
                            value=inorgCompounds.count())

    #inorganic:water molar ratio

    waterMols = DRP.models.CompoundQuantity.objects.filter(reaction=reaction, compound__CSID=937).aggregate(Sum('amount')) 

    if sumInorgAmount > 0 and waterMols > 0:
        num.objects.get_or_create(
                                reaction = reaction,
                                descriptor = descriptorDict['inorgWaterMolRatio'],
                                value = sumInorgAmount/waterMols)

    #inorganic:organic molar ratio

    organicMols = DRP.models.CompoundQuantity.objects.filter(reaction=reaction, role__label='Org'.aggregate(Sum('amount'))
    
    if sumInorgAmount > 0 and organicMols > 0:
        num.objects.get_or_create(
                                reaction=reaction,
                                descriptor=descriptorDict['inorgOrgMolRatio'],
                                value = sumInorgAmount/organicMols)

    #notwater:water molar ratio

    notWaterMols = DRP.models.CompoundQuantity.objects.filter(reaction=reaction).exclude(compound__CSID=937).aggregate(Sum('amount'))
    
    if notWaterMols > 0 and waterMols > 0:
        num.objects.get_or_create(
                                reaction=reaction,
                                descriptor=descriptorDict('notWaterWaterMolRatio'),
                                value = notwaterMols/waterMols 

    for prop in inorgAtomicProperties:
        #  Calculate the inorganic atomic properties, weight them where needed and insert them into the database.
        stem = 'drpInorgAtom' + prop.title().replace('_', '')

        num.objects.get_or_create(
                                reaction=reaction,
                                descriptor=descriptorDict[stem + 'Max'],
                                value=inorgAtomicAggregate(prop, MAX))
    
        num.objects.get_or_create(
                                reaction=reaction,
                                descriptor=descriptorDict[stem + 'Range'],
                                value=inorgAtomicAggregate(prop, RANGE))
    
        num.objects.get_or_create(
                                reaction=reaction,
                                descriptor=descriptorDict[stem + 'Geom'],
                                value=inorgAtomicAggregate(prop, GMEAN))
       
        num.objects.get_or_create(
                                reaction=reaction,
                                descriptor=descriptorDict[stem + 'MolWeightedMax'],
                                value=inorgAtomicAggregate(prop, MAX, WMOL)) 
    
        num.objects.get_or_create(
                                reaction=reaction,
                                descriptor=descriptorDict[stem + 'MolWeightedRange'],
                                value=inorgAtomicAggregate(prop, RANGE, WMOL))
    
        num.objects.get_or_create(
                                reaction=reaction,
                                descriptor=descriptorDict[stem + 'MolWeightedGeom'],
                                value=inorgAtomicAggregate(prop, GMEAN, WMOL))
    
        num.objects.get_or_create(
                                reaction=reaction,
                                descriptor=descriptorDict[stem + 'StoichWeightedMax'],
                                value=inorgAtomicAggregate(prop, MAX, WMOL))
    
        num.objects.get_or_create(
                                reaction=reaction,
                                descriptor=descriptorDict[stem + 'StoichWeightedRange'],
                                value=inorgAtomicAggregate(prop, RANGE, WMOL))
    
        num.objects.get_or_create(
                                reaction=reaction,
                                descriptor=descriptorDict[stem + 'StoichWeightedGeom'],
                                value=inorgAtomicAggregate(prop, GMEAN, WSTOICH))

        # Calculate the elemental molarities
        allCompoundQuantities = CompoundQuantity.objects.filter(reaction=reaction)

        elementNormalisationFactor = sum(sum(quantity.compound.elements[element]['stoichiometry'] * quantity.amount for quantity in allCompoundQuantities) for element in elements)
        # This has been spelled with an s. Someone English has been here...

        for element in elements:
            num.objects.get_or_create(
                                reaction=reaction,
                                descriptor=descriptorDict[element + '_mols'],
                                value=sum(quantity.compound.elements[element]['stoichiometry'] * quantity.amount for quantity in allCompoundQuantities)/elementNormalisationFactor)


