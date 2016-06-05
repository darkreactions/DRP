"""An example molecular descriptor plugin to demonstrate the 'shape' that the API requires."""
# I wanted to name this module rdkit, but then we get name conflicts... lol python
from utils import setup
import rdkit.Chem
import DRP

calculatorSoftware = 'DRP_rdkit'

_descriptorDict = {
    'mw': {'type': 'num', 'name': 'Molecular Weight', 'calculatorSoftware': calculatorSoftware, 'calculatorSoftwareVersion': 0.02, 'maximum': None, 'minimum': 0},
}

descriptorDict = setup(_descriptorDict)

pt = rdkit.Chem.GetPeriodicTable()


def calculate_many(compound_set, verbose=False, whitelist=None):
    for i, compound in enumerate(compound_set):
        if verbose:
            print "{}; Compound {} ({}/{})".format(compound, compound.pk, i + 1, len(compound_set))
        calculate(compound, verbose=verbose, whitelist=whitelist)


def calculate(compound, verbose=False, whitelist=None):
    """Calculate the descriptors from this plugin for a compound."""
    heading = 'mw'
    if whitelist is None or heading in whitelist:
        mw = sum(pt.GetAtomicWeight(pt.GetAtomicNumber(str(element))) * float(compound.elements[element]['stoichiometry']) for element in compound.elements)

        v = DRP.models.NumMolDescriptorValue.objects.update_or_create(defaults={'value': mw}, descriptor=descriptorDict[heading], compound=compound)[0]

        try:
            v.full_clean()
        except ValidationError as e:
            warnings.warn('Value {} for compound {} and descriptor {} failed validation. Value set to None. Validation error message: {}'.format(v.value, v.compound, v.descriptor, e.message))
            v.value = None
            v.save()
