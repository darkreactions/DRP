"""An example molecular descriptor plugin to demonstrate the 'shape' that the API requires."""
from django.conf import settings
from chemspipy import ChemSpider
from utils import setup
import rdkit.Chem
import DRP

_descriptorDict = {
    'mw': {'type': 'num', 'name': 'Molecular Weight', 'calculatorSoftware': 'drp_rdkit', 'calculatorSoftwareVersion': 0, 'maximum': None, 'minimum': 0},
    'fs': {'type': 'ord', 'name': 'Fake size', 'calculatorSoftware': 'example_plugin', 'calculatorSoftwareVersion': 0, 'maximum': 3, 'minimum': 1},
    'N?': {'type': 'bool', 'name': 'Has Nitrogen', 'calculatorSoftware': 'example_plugin', 'calculatorSoftwareVersion': 0},
    'arb': {'type': 'cat', 'name': "Phil's arbitrary descriptor", 'calculatorSoftware': 'example_plugin', 'calculatorSoftwareVersion': 0, 'permittedValues': ('fun', 'dull')}
}
"""A dictionary describing the descriptors available in this module. The key should always be the heading for the descriptor."""

descriptorDict = setup(_descriptorDict)


def fsValueCalc(mw):
    """Calculate an ordinal fake size value."""
    if mw < 50:
        return 1
    elif mw < 100:
        return 2
    else:
        return 3

def calculate_many(compound_set, verbose=False):
    for compound in compound_set:
        calculate(compound)

def arbValCalc(compound):
    """Calculate a completely arbitrary value as an example of a categorical descriptor."""
    if compound.pk % 2 == 0:
        return DRP.models.CategoricalDescriptorPermittedValue.objects.get(value='dull', descriptor=descriptorDict['arb'])
    else:
        return DRP.models.CategoricalDescriptorPermittedValue.objects.get(value='fun', descriptor=descriptorDict['arb'])


def calculate(compound):
    """Calculate the descriptors from this plugin for a compound.

    This should fail silently if a descriptor cannot be calculated for a compound, storing a None value in the
    database as this happens.
    """
    pt = rdkit.Chem.GetPeriodicTable()
    mwValue = DRP.models.NumMolDescriptorValue.objects.get_or_create(descriptor=descriptorDict['mw'], compound=compound)[0]
    mwValue.value = sum(pt.GetAtomicWeight(pt.GetAtomicNumber(str(element))) * float(compound.elements[element]['stoichiometry']) for element in compound.elements)
    mwValue.save()
    fsValue = DRP.models.OrdMolDescriptorValue.objects.get_or_create(compound=compound, descriptor=descriptorDict['fs'])[0]
    fsValue.value = fsValueCalc(mwValue)
    fsValue.save()
    if compound.smiles is None:
        nValue = DRP.models.BoolMolDescriptorValue.objects.get_or_create(compound=compound, descriptor=descriptorDict['N?'])[0]
        nValue.value = None
    else:
        nValue = DRP.models.BoolMolDescriptorValue.objects.get_or_create(compound=compound, descriptor=descriptorDict['N?'])[0]
        nValue.value = ('n' in compound.smiles or 'N' in compound.smiles)
    nValue.save()
    arbValue = DRP.models.CatMolDescriptorValue.objects.get_or_create(compound=compound, descriptor=descriptorDict['arb'])[0]
    arbValue.value = arbValCalc(compound)
    arbValue.save()
