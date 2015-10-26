from scipy.stats import gmean
from utils import setup
import DRP

elements = DRP.chemical_data.elements

#Inorganic descriptors
inorgAtomicProperties = (
    'ionization_energy',
    'electron_affinity',
    'pauling_electronegativity',
    'pearson_electronegativity',
    'hardness',
    'atomic_radius'
)

weightings= (
    ('unw', 'unweighted'),
    ('stoich', 'stoichiometry')
)

inorgElements = {}
for element, info in elements.items():
    if (element == 'Se') or (info['group'] in range(3, 13)) or ((info['group'] > 12) and ((not info['nonmetal']) or info['metalloid'])):
        inorgElements[element] = info 

for prop in inorgAtomicProperties:
    stem = 'drpInorgAtom' + prop.title().replace('_', '') 
    for weighting in weightings:
        _descriptorDict['{}_geom_{}'.format(stem, weighting[0])] = {
            'type':'num',
            'name': 'Geometric mean of {} weighted by {}.'.format(agg[1], prop.replace('_', ' '), weighting[1]),
            'calculatorSoftware':'DRP',
            'calculatorSoftwareVersion':'0.02',
            'maximum':None,
            'minimum':None
            }
    _descriptorDict['{}_max'.format(stem)] = {
        'type': 'num',
        'name': 'Maximal value of {}'.format(prop.replace('_', '')),
        'calculatorSoftware':'DRP',
        'calculatorSoftwareVersion':'0.02',
        'maximum':None,
        'minimum':None
        }
    _descriptorDict['{}_range'.format(stem)] = {
        'type': 'num',
        'name': 'Range of {}'.format(prop.replace('_', '')),
        'calculatorSoftware':'DRP',
        'calculatorSoftwareVersion':'0.02',
        'maximum':None,
        'minimum':None
        }

descriptorDict = setup(_descriptorDict)

def calculate(compound):

    num = DRP.models.NumMolDescriptorValue
    inorgElementNormalisationFactor = sum(info['stoichiometry'] for element, info in compound.elements.items() if element in inorgElements)

    for prop in inorgAtomicProperties:
        num.objects.get_or_create(        
            compound=compound,
            descriptor=descriptorDict['drpInorgAtom{}_geom_unw'.format(prop.title().replace('_', ''))]
            value = gmean(inorgElements[element][prop] for element in compound.elements if element in inorgElements))

        num.objects.get_or_create(
            compound=compound
            descriptor=descriptorDict['drpInorgAtom{}_geom_stoich'.format(prop.title().replace('_', ''))]
            value = gmean(inorgElements[element][prop]*(info['stoichiometry']/inorgElementNormalisationFactor) for element, info in compound.elements if element in inorgElements))
           
        num.objects.get_or_create(
            compound=compound
            descriptor=descriptorDict['drpInorgAtom{}_max'.format(prop.title().replace('_', ''))]
            value = max(inorgElements[element][prop] for element in compound.elements if element in inorgElements))

        num.objects.get_or_create(
            compound=compound
            descriptor=descriptorDict['drpInorgAtom{}_range'.format(prop.title().replace('_', ''))]
            value = max(inorgElements[element][prop] for element in compound.elements if element in inorgElements) - min(inorgElements[element][prop] for element in compound.elements if element in inorgElement))
