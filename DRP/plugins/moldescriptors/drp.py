from scipy.stats import gmean
from utils import setup
import DRP
from DRP import chemical_data
import warnings
from django.core.exceptions import ValidationError

calculatorSoftware = 'DRP'

create_threshold = 5000  # number of values to create at a time. Should probably be <= 5000

elements = chemical_data.elements

# Inorganic descriptors
inorgAtomicProperties = (
    'ionization_energy',
    'electron_affinity',
    'pauling_electronegativity',
    'pearson_electronegativity',
    'hardness',
    'atomic_radius'
)

weightings = (
    ('unw', 'unweighted'),
    ('stoich', 'stoichiometry')
)

inorgElements = {}
for element, info in elements.items():
    if (element == 'Se') or (info['group'] in range(3, 13)) or ((info['group'] > 12) and ((not info['nonmetal']) or info['metalloid'])):
        inorgElements[element] = info

_descriptorDict = {}

for prop in inorgAtomicProperties:
    stem = 'drpInorgAtom' + prop.title().replace('_', '')
    for weighting in weightings:
        _descriptorDict['{}_geom_{}'.format(stem, weighting[0])] = {
            'type': 'num',
            'name': 'Geometric mean of {} weighted by {}.'.format(prop.replace('_', ' '), weighting[1]),
            'calculatorSoftware': calculatorSoftware,
            'calculatorSoftwareVersion': '0_02',
            'maximum': None,
            'minimum': 0,
        }
    _descriptorDict['{}_max'.format(stem)] = {
        'type': 'num',
        'name': 'Maximal value of {}'.format(prop.replace('_', '')),
        'calculatorSoftware': calculatorSoftware,
        'calculatorSoftwareVersion': '0_02',
        'maximum': None,
        'minimum': None
    }
    _descriptorDict['{}_range'.format(stem)] = {
        'type': 'num',
        'name': 'Range of {}'.format(prop.replace('_', '')),
        'calculatorSoftware': calculatorSoftware,
        'calculatorSoftwareVersion': '0_02',
        'maximum': None,
        'minimum': 0,
    }

for group_num in range(1, 19):
    _descriptorDict['drpInorgAtom_boolean_group_{}'.format(group_num)] = {
        'type': 'bool',
        'name': 'Presence of inorganic elements in group {}'.format(group_num),
        'calculatorSoftware': calculatorSoftware,
        'calculatorSoftwareVersion': '1_5',
    }

for period_num in range(1, 8):
    _descriptorDict['drpInorgAtom_boolean_period_{}'.format(period_num)] = {
        'type': 'bool',
        'name': 'Presence of inorganic elements in period {}'.format(period_num),
        'calculatorSoftware': calculatorSoftware,
        'calculatorSoftwareVersion': '1_5',
    }

for valence_num in range(0, 8):
    _descriptorDict['drpInorgAtom_boolean_valence_{}'.format(valence_num)] = {
        'type': 'bool',
        'name': 'Presence of inorganic elements with valence {}'.format(valence_num),
        'calculatorSoftware': calculatorSoftware,
        'calculatorSoftwareVersion': '1_5',
    }

for group_num in range(1, 19):
    _descriptorDict['boolean_group_{}'.format(group_num)] = {
        'type': 'bool',
        'name': 'Presence of elements in group {}'.format(group_num),
        'calculatorSoftware': calculatorSoftware,
        'calculatorSoftwareVersion': '1_5',
    }

for period_num in range(1, 8):
    _descriptorDict['boolean_period_{}'.format(period_num)] = {
        'type': 'bool',
        'name': 'Presence of elements in period {}'.format(period_num),
        'calculatorSoftware': calculatorSoftware,
        'calculatorSoftwareVersion': '1_5',
    }

descriptorDict = setup(_descriptorDict)


def delete_descriptors(compound_set, whitelist=None):
    if whitelist is None:
        descs = descriptorDict.values()
    else:
        descs = [descriptorDict[k] for k in descriptorDict.keys() if k in whitelist]
    DRP.models.NumMolDescriptorValue.objects.filter(descriptor__in=[desc for desc in descs if isinstance(desc, DRP.models.NumMolDescriptor)], compound__in=compound_set).delete(recalculate_reactions=False)
    DRP.models.BoolMolDescriptorValue.objects.filter(descriptor__in=[desc for desc in descs if isinstance(desc, DRP.models.BoolMolDescriptor)], compound__in=compound_set).delete(recalculate_reactions=False)


def calculate_many(compound_set, verbose=False, whitelist=None):
    delete_descriptors(compound_set, whitelist=whitelist)
    num_vals_to_create = []
    bool_vals_to_create = []
    for i, compound in enumerate(compound_set):
        if verbose:
            print "{}; Compound {} ({}/{})".format(compound, compound.pk, i + 1, len(compound_set))
        num_vals_to_create, bool_vals_to_create = _calculate(compound, num_vals_to_create=num_vals_to_create, bool_vals_to_create=bool_vals_to_create)
        if len(num_vals_to_create) > create_threshold:
            if verbose:
                print 'Creating {} numeric values'.format(len(num_vals_to_create))
                DRP.models.NumMolDescriptorValue.objects.bulk_create(num_vals_to_create)
                num_vals_to_create = []
        if len(bool_vals_to_create) > create_threshold:
            if verbose:
                print 'Creating {} boolean values'.format(len(bool_vals_to_create))
            DRP.models.BoolMolDescriptorValue.objects.bulk_create(bool_vals_to_create)
            bool_vals_to_create = []

    if verbose:
        print 'Creating {} numeric values'.format(len(num_vals_to_create))
    DRP.models.NumMolDescriptorValue.objects.bulk_create(num_vals_to_create)
    if verbose:
        print 'Creating {} boolean values'.format(len(bool_vals_to_create))
    DRP.models.BoolMolDescriptorValue.objects.bulk_create(bool_vals_to_create)


def calculate(compound):
    delete_descriptors([compound])
    num_vals_to_create, bool_vals_to_create = _calculate(compound)
    DRP.models.NumMolDescriptorValue.objects.bulk_create(num_vals_to_create)
    if verbose:
        print 'Creating {} numeric values'.format(len(num_vals_to_create))
    DRP.models.NumMolDescriptorValue.objects.bulk_create(num_vals_to_create)
    if verbose:
        print 'Creating {} boolean values'.format(len(bool_vals_to_create))
    DRP.models.BoolMolDescriptorValue.objects.bulk_create(bool_vals_to_create)


def _calculate(compound, verbose=False, whitelist=None, num_vals_to_create=[], bool_vals_to_create=[]):

    num = DRP.models.NumMolDescriptorValue
    boolVal = DRP.models.BoolMolDescriptorValue

    if any(element in inorgElements for element in compound.elements):
        inorgElementNormalisationFactor = sum(info['stoichiometry'] for element, info in compound.elements.items() if element in inorgElements)
        for prop in inorgAtomicProperties:
            heading = 'drpInorgAtom{}_geom_unw'.format(prop.title().replace('_', ''))
            if whitelist is None or heading in whitelist:
                # zero is what scipy does natively. This is just to avoid warnings that are fine so they don't drown out the real ones
                if any([(inorgElements[element][prop] == 0) for element, info in compound.elements.items() if element in inorgElements]):
                    val = 0
                elif any([(inorgElements[element][prop] < 0) for element, info in compound.elements.items() if element in inorgElements]):
                    raise ValueError('Cannot take geometric mean of negative values. This descriptor ({}) should not use a geometric mean.'.format(descriptorDict['drpInorgAtom{}_geom_unw'.format(prop.title().replace('_', ''))]))
                else:
                    val = gmean([inorgElements[element][prop] for element in compound.elements if element in inorgElements])
                n = num(
                    compound=compound,
                    descriptor=descriptorDict[heading],
                    value=val
                )
                try:
                    n.full_clean()
                except ValidationError as e:
                    warnings.warn('Value {} for compound {} and descriptor {} failed validation. Value set to none. Validation error message: {}'.format(n.value, n.compound, n.descriptor, e.message))
                    n.value = None
                num_vals_to_create.append(n)

            heading = 'drpInorgAtom{}_geom_stoich'.format(prop.title().replace('_', ''))
            if whitelist is None or heading in whitelist:
                if any([(inorgElements[element][prop] * float(info['stoichiometry'] / inorgElementNormalisationFactor) == 0) for element, info in compound.elements.items() if element in inorgElements]):
                    val = 0
                elif any([(inorgElements[element][prop] * float(info['stoichiometry'] / inorgElementNormalisationFactor) < 0) for element, info in compound.elements.items() if element in inorgElements]):
                    raise ValueError('Cannot take geometric mean of negative values. This descriptor ({}) should not use a geometric mean.'.format(descriptorDict['drpInorgAtom{}_geom_stoich'.format(prop.title().replace('_', ''))]))
                else:
                    try:
                        val = gmean([inorgElements[element][prop] * float(info['stoichiometry'] / inorgElementNormalisationFactor) for element, info in compound.elements.items() if element in inorgElements])
                    except:
                        print [inorgElements[element][prop] * float((info['stoichiometry'] / inorgElementNormalisationFactor)) for element, info in compound.elements.items() if element in inorgElements]
                        print [inorgElements[element][prop] for element, info in compound.elements.items() if element in inorgElements]
                        print [float((info['stoichiometry'] / inorgElementNormalisationFactor)) for element, info in compound.elements.items() if element in inorgElements]
                        print [info['stoichiometry'] / inorgElementNormalisationFactor for element, info in compound.elements.items() if element in inorgElements]
                        print [info['stoichiometry'] for element, info in compound.elements.items() if element in inorgElements]
                        print [inorgElementNormalisationFactor for element, info in compound.elements.items() if element in inorgElements]
                        raise
                n = num(
                    compound=compound,
                    descriptor=descriptorDict[heading],
                    value=val,
                )
                try:
                    n.full_clean()
                except ValidationError as e:
                    warnings.warn('Value {} for compound {} and descriptor {} failed validation. Value set to none. Validation error message: {}'.format(n.value, n.compound, n.descriptor, e.message))
                    n.value = None
                num_vals_to_create.append(n)

            heading = 'drpInorgAtom{}_max'.format(prop.title().replace('_', ''))
            if whitelist is None or heading in whitelist:
                val = max(inorgElements[element][prop] for element in compound.elements if element in inorgElements)
                n = num(
                    compound=compound,
                    descriptor=descriptorDict[heading],
                    value=val
                )
                try:
                    n.full_clean()
                except ValidationError as e:
                    warnings.warn('Value {} for compound {} and descriptor {} failed validation. Value set to none. Validation error message: {}'.format(n.value, n.compound, n.descriptor, e.message))
                    n.value = None
                num_vals_to_create.append(n)

            heading = 'drpInorgAtom{}_range'.format(prop.title().replace('_', ''))
            if whitelist is None or heading in whitelist:
                val = max(inorgElements[element][prop] for element in compound.elements if element in inorgElements) - min(inorgElements[element][prop] for element in compound.elements if element in inorgElements)
                n = num(
                    compound=compound,
                    descriptor=descriptorDict[heading],
                    value=val,
                )
                try:
                    n.full_clean()
                except ValidationError as e:
                    warnings.warn('Value {} for compound {} and descriptor {} failed validation. Value set to none. Validation error message: {}'.format(n.value, n.compound, n.descriptor, e.message))
                    n.value = None
                num_vals_to_create.append(n)

    # inorganic elements only
    for group_num in range(1, 19):
        heading = 'drpInorgAtom_boolean_group_{}'.format(group_num)
        if whitelist is None or heading in whitelist:
            bool_vals_to_create.append(boolVal(
                compound=compound,
                descriptor=descriptorDict[heading],
                value=(any(elements[element]['group'] == group_num for element in compound.elements.keys() if element in inorgElements))
            ))
    for period_num in range(1, 8):
        heading = 'drpInorgAtom_boolean_period_{}'.format(period_num)
        if whitelist is None or heading in whitelist:
            bool_vals_to_create.append(boolVal(
                compound=compound,
                descriptor=descriptorDict[heading],
                value=(any(elements[element]['period'] == period_num for element in compound.elements.keys() if element in inorgElements))
            ))
    for valence_num in range(0, 8):
        heading = 'drpInorgAtom_boolean_valence_{}'.format(valence_num)
        if whitelist is None or heading in whitelist:
            bool_vals_to_create.append(boolVal(
                compound=compound,
                descriptor=descriptorDict[heading],
                value=(any(elements[element]['valence'] == valence_num for element in compound.elements.keys() if element in inorgElements))
            ))
    # all elements
    for group_num in range(1, 19):
        heading = 'boolean_group_{}'.format(group_num)
        if whitelist is None or heading in whitelist:
            bool_vals_to_create.append(boolVal(
                compound=compound,
                descriptor=descriptorDict[heading],
                value=(any(elements[element]['group'] == group_num for element in compound.elements.keys()))
            ))
    for period_num in range(1, 8):
        heading = 'boolean_period_{}'.format(period_num)
        if whitelist is None or heading in whitelist:
            bool_vals_to_create.append(boolVal(
                compound=compound,
                descriptor=descriptorDict[heading],
                value=(any(elements[element]['period'] == period_num for element in compound.elements.keys()))
            ))

    return num_vals_to_create, bool_vals_to_create
