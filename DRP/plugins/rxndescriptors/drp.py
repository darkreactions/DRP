"""Basic reaction descriptors calculation module"""

from django.conf import settings
from utils import setup
import DRP 

atoms = DRP.chemical_data.atoms

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
}

descriptorDict = setup(_descriptorDict)

def calculate(reaction):
    num = DRP.models.NumRxnDescriptorValue
    num.objects.get_or_create(
                            reaction=reaction,
                            descriptor=descriptorDict['numberInorganic'],
                            value=DRP.models.CompoundRole.objects.filter(reaction=reaction, role__label='Inorg').count())

    num.objects.get_or_create(
                            reaction=reaction,
                            descriptor=descriptorDict['inorgAtomIonizationMax'],
                            value=max(max(atoms[element]['ionization_energy'] for element in compound.elements) for compound in reaction.compounds))

    num.objects.get_or_create(
                            reaction=reaction,
                            descriptor=descriptorDict['inorgAtomIonizationMin'],
                            value=min(min(atoms[element]['ionization_energy'] for element in compound.elements) for compound in reaction.compounds))
