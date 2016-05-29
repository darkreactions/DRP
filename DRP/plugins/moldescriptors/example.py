"""An example molecular descriptor plugin to demonstrate the 'shape' that the API requires."""
from utils import setup
import DRP

_descriptorDict = {
    'fs': {'type': 'ord', 'name': 'Fake size', 'calculatorSoftware': 'example_plugin', 'calculatorSoftwareVersion': 1.5, 'maximum': 3, 'minimum': 1},
    'N?': {'type': 'bool', 'name': 'Has Nitrogen', 'calculatorSoftware': 'example_plugin', 'calculatorSoftwareVersion': 0},
    'arb': {'type': 'cat', 'name': "Phil's arbitrary descriptor", 'calculatorSoftware': 'example_plugin', 'calculatorSoftwareVersion': 0, 'permittedValues': ('fun', 'dull')}
}
"""A dictionary describing the descriptors available in this module. The key should always be the heading for the descriptor."""

descriptorDict = setup(_descriptorDict)


def fsValueCalc(num):
    """Calculate an ordinal fake size value."""
    if num < 10:
        return 1
    elif num < 20:
        return 2
    else:
        return 3


def arbValCalc(compound):
    """Calculate a completely arbitrary value as an example of a categorical descriptor."""
    if compound.pk % 2 == 0:
        return DRP.models.CategoricalDescriptorPermittedValue.objects.get(value='dull', descriptor=descriptorDict['arb'])
    else:
        return DRP.models.CategoricalDescriptorPermittedValue.objects.get(value='fun', descriptor=descriptorDict['arb'])


def calculate_many(compound_set, verbose=False, whitelist=None):
    for i, compound in enumerate(compound_set):
        if verbose:
            print "{}; Compound {} ({}/{})".format(compound, compound.pk, i + 1, len(compound_set))
        calculate(compound, verbose=verbose, whitelist=whitelist)


def calculate(compound, verbose=False, whitelist=None):
    """Calculate the descriptors from this plugin for a compound.

    This should fail silently if a descriptor cannot be calculated for a compound, storing a None value in the
    database as this happens.
    """

    if compound.smiles is None:
        fsValue = None
        nValue = None
    else:
        nValue = ('n' in compound.smiles or 'N' in compound.smiles)
        fsValue = fsValueCalc(len(compound.smiles))

    arbValue = arbValCalc(compound)
    heading = 'fs'
    if whitelist is None or heading in whitelist:
        v = DRP.models.OrdMolDescriptorValue.objects.update_or_create(defaults={'value': fsValue}, compound=compound, descriptor=descriptorDict[heading])[0]
        try:
            v.full_clean()
        except ValidationError as e:
            warnings.warn('Value {} for compound {} and descriptor {} failed validation. Value set to None. Validation error message: {}'.format(v.value, v.compound, v.descriptor, e.message))
            v.value = None
            v.save()
            
    heading = 'N?'
    if whitelist is None or heading in whitelist:
        v = DRP.models.BoolMolDescriptorValue.objects.update_or_create(defaults={'value': nValue}, compound=compound, descriptor=descriptorDict[heading])[0]
        try:
            v.full_clean()
        except ValidationError as e:
            warnings.warn('Value {} for compound {} and descriptor {} failed validation. Value set to None. Validation error message: {}'.format(v.value, v.compound, v.descriptor, e.message))
            v.value = None
            v.save()

    heading = 'arb'
    if whitelist is None or heading in whitelist:
        v = DRP.models.CatMolDescriptorValue.objects.update_or_create(defaults={'value': arbValue}, compound=compound, descriptor=descriptorDict[heading])[0]
        try:
            v.full_clean()
        except ValidationError as e:
            warnings.warn('Value {} for compound {} and descriptor {} failed validation. Value set to None. Validation error message: {}'.format(v.value, v.compound, v.descriptor, e.message))
            v.value = None
            v.save()
    
