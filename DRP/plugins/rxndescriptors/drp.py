"""Basic reaction descriptors calculation module"""
from itertools import chain
from numpy import mean, average as wmean
from scipy.stats import gmean
from django.conf import settings
from utils import setup
import DRP 

atoms = DRP.chemical_data.elements

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
    ('Min', 'Minimum'),
    ('Mean', 'Mean Average'),
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

#Set up the actual descriptor dictionary.
descriptorDict = setup(_descriptorDict)

for element, info in elements.items():
    if (element == 'Se') or (info['group'] in range(3, 13)) or ((info['group'] > 12) and ((not info['nonmetal']) or info['metalloid'])):
        inorgElements[element] = info 


def calculate(reaction):
    """Calculate the descriptors for this plugin."""
    num = DRP.models.NumRxnDescriptorValue
    inorgCompoundQuantities = DRP.models.CompoundQuantity.objects.filter(reaction=reaction, role__label='Inorg')  # These two lines get the inorganic compounds for htis reaction.
    sumInorgAmount = sum(quantity.amount for quantity in inorgCompoundQuantities)


    MAX = 1
    MIN = 2
    MEAN = 3
    GMEAN = 4
    WMOL = 1
    WSTOICH = 2

    def inorgAtomicAggregate(propertyLabel, function, weighting=None):
        """Figure out the aggregation for inorganic atom-based properties. This is done a lot so this saves on code repetition."""
        if weighting is None:
            if function == MAX:
                return max(max(inorgElements[element][propertyLabel] for element in quantity.compound.elements) for quantity in inorgCompoundQuantities)
            elif function == MIN:
                return min(min(inorgElements[element][propertyLabel] for element in quantity.compound.elements) for quantity in inorgCompoundQuantities)
            elif function == MEAN:
                return mean(chain(inorgElements[element][propertyLabel] for element in quantity.compound.elements) for quantity in inorgCompoundQuantities)
            elif function == GMEAN:
                return gmean(chain(inorgElements[element][propertyLabel] for element in quantity.compound.elements) for quantity in inorgCompoundQuantities)
            else:
                raise TypeError('Unrecognised function selection constant')
        elif weighting == WMOL:
            if function == MAX:
                return max(max(inorgElements[element][propertyLabel]*(quantity.amount/sumInorgAmount)*info['stoichiometry'] for element, info in quantity.compound.elements)
                                            for quantity in inorgCompoundQuantities)
            elif function == MIN:
                return min(min(inorgElements[element][propertyLabel]*(quantity.amount/sumInorgAmount)*info['stoichiometry'] for element, info in quantity.compound.elements)
                                            for quantity in inorgCompoundQuantities)
            elif function == MEAN:
                return mean(chain(inorgElements[element][propertyLabel]*(quantity.amount/sumInorgAmount)*info['stoichiometry'] for element, info in quantity.compound.elements)
                                            for quantity in inorgCompoundQuantities)
            elif function == GMEAN:
                return gmean(chain(inorgElements[element][propertyLabel]*(quantity.amount/sumInorgAmount)*info['stoichiometry'] for element, info in quantity.compound.elements)
                                            for quantity in inorgCompoundQuantities))
            else:
                raise TypeError('Unrecognised function selection constant')
        elif weighting == WSTOICH: 
            values = (wmean((inorgElements[element][propertyLabel] for element in quantity.compound.elements), 
                                weights=(info['stoichiometry'] for element, info in quantity.compound.elements.items()))
                                            for quantity in inorgCompoundQuantities)
            if function == MAX:
                return max(values)
            elif function == MAX:
                return min(values)
            elif function == MEAN:
                return mean(values)
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

    for prop in inorgAtomicProperties:
        #  Calculate the inorganic atomic properties, weight them where needed and insert them into the database.
        stem = 'drpInorgAtom' + prop.title().replace('_', '')

        num.objects.get_or_create(
                                reaction=reaction,
                                descriptor=descriptorDict[stem + 'Max'],
                                value=inorgAtomicAggregate(prop, MAX))
    
        num.objects.get_or_create(
                                reaction=reaction,
                                descriptor=descriptorDict[stem + 'Min'],
                                value=inorgAtomicAggregate(prop, MIN))
    
        num.objects.get_or_create(
                                reaction=reaction,
                                descriptor=descriptorDict[stem + 'Mean'],
                                value=inorgAtomicAggregate(prop, MEAN))
    
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
                                descriptor=descriptorDict[stem + 'MolWeightedMax'],
                                value=inorgAtomicAggregate(prop, MIN, WMOL))
    
        num.objects.get_or_create(
                                reaction=reaction,
                                descriptor=descriptorDict[stem + 'MolWeightedMean'],
                                value=inorgAtomicAggregate(prop, MEAN, WMOL))
    
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
                                descriptor=descriptorDict[stem + 'StoichWeightedMin'],
                                value=inorgAtomicAggregate(prop, MIN, WMOL))
    
        num.objects.get_or_create(
                                reaction=reaction,
                                descriptor=descriptorDict[stem + 'StoichWeightedMean'],
                                value=inorgAtomicAggregate(prop, MEAN, WSTOICH))
    
        num.objects.get_or_create(
                                reaction=reaction,
                                descriptor=descriptorDict[stem + 'StoichWeightedGeom'],
                                value=inorgAtomicAggregate(prop, GMEAN, WSTOICH))
