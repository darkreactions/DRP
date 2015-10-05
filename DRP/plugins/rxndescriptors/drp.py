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
        }
    'inorgAtomIonizationMax':
        {
            'type': 'num',
            'name': 'Maximum atomic ionization energy among inorganic species',
            'calculatorSoftware': 'DRP',
            'calculatorSoftwareversion':'0.02',
            'maximum': None,
            'minimum': None
        }
    'inorgAtomIonizationMin':
        {
            'type': 'num',
            'name': 'Minimum atomic ionization energy among inorganic species',
            'calculatorSoftware': 'DRP',
            'calculatorSoftwareversion':'0.02',
            'maximum': None,
            'minimum': None
        }
    'inorgAtomIonizationMean':
        {
            'type': 'num',
            'name': 'Mean atomic ionization energy among inorganic species',
            'calculatorSoftware': 'DRP',
            'calculatorSoftwareversion':'0.02',
            'maximum': None,
            'minimum': None
        }
    'inorgAtomIonizationGeom'
        {
            'type': 'num',
            'name': 'Geometric average atomic ionization energy among inorganic species',
            'calculatorSoftware': 'DRP',
            'calculatorSoftwareversion':'0.02',
            'maximum': None,
            'minimum': None
        }
    'inorgAtomIonizationMaxWeighted'
        {
            'type': 'num',
            'name': 'Maximum atomic ionization energy among inorganic species, stoichiometrically and molar weighted.',
            'calculatorSoftware': 'DRP',
            'calculatorSoftwareversion':'0.02',
            'maximum': None,
            'minimum': None
        }
}

descriptorDict = setup(_descriptorDict)

InorgElements = {}

for element, info in elements.items():
    if (element == 'Se') or (info['group'] in range(3, 13)) or ((info['group'] > 12) and ((not info['nonmetal']) or info['metalloid'])):
        inorgElements[element] = info 


def calculate(reaction):
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
                return mean(chain(inorgElements[element]['ionization_energy']*(quantity.amount/sumInorgAmount)*info['stoichiometry'] for element, info in quantity.compound.elements)
                                            for quantity in inorgCompoundQuantities)
            elif function == GMEAN:
                return gmean(chain(inorgElements[element]['ionization_energy']*(quantity.amount/sumInorgAmount)*info['stoichiometry'] for element, info in quantity.compound.elements)
                                            for quantity in inorgCompoundQuantities))
            else:
                raise TypeError('Unrecognised function selection constant')
        elif weighting == WSTOICH: 
            values = (wmean((inorgElements[element]['ionization_energy'] for element in quantity.compound.elements), 
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

    # Maximum ionisation potential across all relevant inorganic atom types in the reaction.
    num.objects.get_or_create(
                            reaction=reaction,
                            descriptor=descriptorDict['drpInorgMetalAtomIonizationMax'],
                            value=inorgAtomicAggregate('ionization_energy', MAX)

    # Minimum ionisation potential across all relevant inorganic atom types in the reaction.
    num.objects.get_or_create(
                            reaction=reaction,
                            descriptor=descriptorDict['drpInorgMetalAtomIonizationMin'],
                            value=inorgAtomicAggregate('ionization_energy', MIN)

    # Mean ionisation potential across all relevant inorganic atom types in the reaction
    num.objects.get_or_create(
                            reaction=reaction,
                            descriptor=descriptorDict['drpInorgMetalAtomIonizationMean'],
                            value=inorgAtomicAggregate('ionization_energy', MEAN)

    # Geometric Average of the ionisation potential across all relevant inorganic atom types in the reaction
    num.objects.get_or_create(
                            reaction=reaction,
                            descriptor=descriptorDict['drpInorgMetalAtomIonizationGeom'],
                            value=inorgAtomicAggregate('ionization_energy', GMEAN)
   
    # Maximum molarity-weighted value of ionisation potential across all relevant inorganic atoms in the reaction
    num.objects.get_or_create(
                            reaction=reaction,
                            descriptor=descriptorDict['drpInorgMetalAtomIonizationMaxMolWeighted'],
                            value=inorgAtomicAggregate('ionization_energy', MAX, WMOL) 

    # Minimum molarity-weighted value of ionisation potential across all relevant inorganic atoms in the reaction
    num.objects.get_or_create(
                            reaction=reaction,
                            descriptor=descriptorDict['drpInorgMetalAtomIonizationMinMolWeighted'],
                            value=inorgAtomicAggregate('ionization_energy', MIN, WMOL)

    # Average molarity-weighted value of ionisation potential across all relevant inorganic atoms in the reaction
    num.objects.get_or_create(
                            reaction=reaction,
                            descriptor=descriptorDict['drpInorgMetalAtomIonizationMeanMolWeighted'],
                            value=inorgAtomicAggregate('ionization_energy', MEAN, WMOL)

    # Geometric Average molarity-weighted value of ionisation potential across all relevant inorganic atoms in the reaction
    num.objects.get_or_create(
                            reaction=reaction,
                            descriptor=descriptorDict['drpInorgMetalAtomIonizationGeomMolWeighted'],
                            value=inorgAtomicAggregate('ionization_energy', GMEAN, WMOL)

    # Maximum (stoichiometrically weighted average value of ionisation potential across all relevant inorganic atoms in each species) across all species in the reaction
    num.objects.get_or_create(
                            reaction=reaction,
                            descriptor=descriptorDict['drpInorgMetalAtomIonizationMaxStoichWeighted'],
                            value=inorgAtomicAggregate('ionization_energy', MAX, WMOL)
    # Minimum (stoichiometrically weighted average value of ionisation potential across all relevant inorganic atoms in each species) across all species in the reaction
    num.objects.get_or_create(
                            reaction=reaction,
                            descriptor=descriptorDict['drpInorgMetalAtomIonizationMinStoichWeighted'],
                            value=inorgAtomicAggregate('ionization_energy', MIN, WMOL)

    # Mean (stoichiometrically weighted value of ionisation potential across all relevant inorganic atoms in each species) across all species in the reaction
    num.objects.get_or_create(
                            reaction=reaction,
                            descriptor=descriptorDict['drpInorgMetalAtomIonizationMeanStoichWeighted'],
                            value=inorgAtomicAggregate('ionization_energy', MEAN, WSTOICH)

    # Geometric Mean (stoichiometrically weighted value of ionisation potential across all relevant inorganic atoms in each species) across all species in the reaction
    num.objects.get_or_create(
                            reaction=reaction,
                            descriptor=descriptorDict['drpInorgMetalAtomIonizationGeomStoichWeighted'],
                            value=inorgAtomicAggregate('ionization_energy', GMEAN, WSTOICH)
