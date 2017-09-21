"""An example molecular descriptor plugin to demonstrate the 'shape' that the API requires."""
# I wanted to name this module rdkit, but then we get name conflicts...
# lol python
from .utils import setup
import DRP
from DRP import chemical_data
<<<<<<< HEAD
import rdkit.Chem
from rdkit.Chem import Descriptors
=======

# importing from rdkit version 2.2016.03.4

import rdkit.Chem
from rdkit.Chem import rdMolDescriptors
>>>>>>> e08a9d8bcd64b253b8f31062a7cf280d17bb3a0e
import logging
from django.core.exceptions import ValidationError
logger = logging.getLogger('DRP')
tracer = logging.getLogger('DRP.tracer')

calculatorSoftware = 'DRP_rdkit'

inorgElements = {}
for element, info in chemical_data.elements.items():
    if info['group'] is not None:
        if (element == 'Se') or (info['group'] in range(3, 13)) or ((info['group'] > 12) and ((not info['nonmetal']) or info['metalloid'])):
            inorgElements[element] = info

weightings = (
    ('unw', 'unweighted'),
    ('stoich', 'stoichiometry')
)

_descriptorDict = {
    'mw': {'type': 'num', 'name': 'Molecular Weight', 'calculatorSoftware': calculatorSoftware, 'calculatorSoftwareVersion': '0_02', 'maximum': None, 'minimum': 0},
    'rbc': {'type': 'num', 'name': 'Rotatable bond count', 'calculatorSoftware': calculatorSoftware, 'calculatorSoftwareVersion': '1_9', 'maximum': None, 'minimum': 0},
    'Chi0v': {'type': 'num', 'name': 'Zero Order Molecular Valence Connectivity Index', 'calculatorSoftware': calculatorSoftware, 'calculatorSoftwareVersion': '1_9', 'maximum': None, 'minimum': 0}
}

elementPropertyKeys = []
for i in range(0, 9):
    for element in inorgElements:
        key = '{}@{}'.format(element, i)
        elementPropertyKeys.append(key)
        _descriptorDict[key] = {
            'type': 'bool',
            'name': 'Presence of element {} at ox state {}'.format(element, i),
            'calculatorSoftware': calculatorSoftware,
            'calculatorSoftwareVersion': '1_9',
        }
tracer.debug('Setting up lazydescdict')
descriptorDict = setup(_descriptorDict)

pt = rdkit.Chem.GetPeriodicTable()


def calculate_many(compound_set, verbose=False, whitelist=None):
    """Calculate in bulk."""
    tracer.debug('Calculating Many')
    whitelist = descriptorDict.keys() if whitelist is None else whitelist
    DRP.models.NumMolDescriptorValue.objects.filter(
        compound__in=compound_set,
        descriptor__in=(descriptor for key, descriptor in descriptorDict.items() if (key in whitelist and key in ('mw', 'rbc', 'Chi0v')))).delete()
    DRP.models.BoolMolDescriptorValue.objects.filter(
        compound__in=compound_set,
        descriptor__in=[descriptor for key, descriptor in descriptorDict.items() if (key in elementPropertyKeys)]).delete()
    resnums = []
    resbools = []
    for i, compound in enumerate(compound_set):
        if verbose:
            logger.info("{}; Compound {} ({}/{})".format(compound,
                                                         compound.pk, i + 1, len(compound_set)))
        nums, bools = calculate(compound, verbose=verbose, whitelist=whitelist)
        resnums += nums
        resbools += bools
    for n in resnums:
        n.save()
#    DRP.models.NumMolDescriptorValue.objects.bulk_create(resnums)
    DRP.models.BoolMolDescriptorValue.objects.bulk_create(resbools)


def validateNumeric(v):
    """Validate a numeric descriptor without breaking the calculation script."""
    try:
        v.clean()
    except ValidationError as e:
        logger.warning('Value {} for compound {} and descriptor {} failed validation. Value set to None. Validation error message: {}'.format(
            v.value, v.compound, v.descriptor, e))
        v.value = None


def recurseSumCharge(atom, already_seen_ids=None):
    """Calculate the formal charge for atoms in an inorganic compound arising from an ingested smiles."""
    if already_seen_ids is None:
        already_seen_ids = set()
    already_seen_ids.add(atom.GetIdx())
    result = atom.GetFormalCharge()
    for neighbor in atom.GetNeighbors():
        if neighbor.GetIdx() not in already_seen_ids:
            result += recurseSumCharge(neighbor, already_seen_ids)
    return result


def calculate(compound, verbose=False, whitelist=None):
    """Calculate the descriptors from this plugin for a compound."""
    heading = 'mw'
    nums = []
    bools = []
    tracer.debug('Calculate function.')
    if whitelist is None or heading in whitelist:
        mw = sum(pt.GetAtomicWeight(pt.GetAtomicNumber(str(element))) * float(
            compound.elements[element]['stoichiometry']) for element in compound.elements)

        v = DRP.models.NumMolDescriptorValue(
            value=mw, descriptor=descriptorDict[heading], compound=compound)
        validateNumeric(v)
        nums.append(v)
    mol = rdkit.Chem.MolFromSmiles(compound.smiles)
    if mol is None:
        logger.warning(
            'Compound {} has no smiles. Skipping calculations for rdkit molecular descriptors.')
    else:
        if whitelist is None or 'rbc' in whitelist:
            rbc = Descriptors.NumRotatableBonds(mol)
            v = DRP.models.NumMolDescriptorValue(
                value=rbc, descriptor=descriptorDict['rbc'], compound=compound)
            validateNumeric(v)
            nums.append(v)
        if whitelist is None or 'Chi0v' in whitelist:
            chi0v = Descriptors.Chi0v(mol)
            v = DRP.models.NumMolDescriptorValue(value=chi0v, descriptor=descriptorDict[
                                                 'Chi0v'], compound=compound)
            validateNumeric(v)
            nums.append(v)
        for element in inorgElements.keys():
            oxStates = []
            for atom in mol.GetAtoms():  # weird capitalisation is weird, but correct.
                if atom.GetSymbol() == element:
                    oxStates.append(recurseSumCharge(
                        atom) + atom.GetTotalValence())
            for ox in range(0, 9):
                oxString = '{}@{}'.format(element, ox)
                if (whitelist is None) or (oxString in whitelist):
                    v = DRP.models.BoolMolDescriptorValue(
                        value=ox in oxStates, descriptor=descriptorDict[oxString], compound=compound)
                    validateNumeric(v)
                    bools.append(v)
    tracer.debug("here are nums: {}".format(str(nums)))
    tracer.debug("here are bools: {}".format(str(bools)))
    return nums, bools
