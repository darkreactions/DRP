"""A descriptor calculator to wrap and use calls to ChemAxon.

Requires cxcalc (part of JChem) to be installed and licensed.

"""
import DRP
from utils import setup
from django.conf import settings
from django.core.exceptions import ValidationError
from collections import OrderedDict
from subprocess import Popen, PIPE
from itertools import chain
import warnings

# The version of ChemAxon currently in use.
# Should be matched by a entry in the dictionary in the settings file
# This version is reflected in, but not necessarily the same as, the
# descriptor versions
CHEMAXON_VERSION = '16.5'

calculatorSoftware = 'ChemAxon_cxcalc'
# number of values to create at a time. Should probably be <= 5000
create_threshold = 5000


# The descriptor versions correspond to either the first ChemAxon version in which they were used
# by the DRP or the oldest ChemAxon version in which a change was made to the methodology for
# that descriptor's calculation (does not include minor bugfixes)

_descriptorDict = {
    'refractivity': {
        'type': 'num',
        'name': 'Refractivity',
        'calculatorSoftware': calculatorSoftware,
        'calculatorSoftwareVersion': '15_6',
        'minimum': 0,
    },
    'maximalprojectionarea': {
        'type': 'num',
        'name': 'Maximal Projection Area',
        'calculatorSoftware': calculatorSoftware,
        'calculatorSoftwareVersion': '15_6',
        'minimum': 0,
    },
    'maximalprojectionradius': {
        'type': 'num',
        'name': 'Maximal Projection Radius',
        'calculatorSoftware': calculatorSoftware,
        'calculatorSoftwareVersion': '15_6',
        'minimum': 0,
    },
    'maximalprojectionsize': {
        'type': 'num',
        'name': 'Maximal Projection Size',
        'calculatorSoftware': calculatorSoftware,
        'calculatorSoftwareVersion': '15_6',
        'minimum': 0,
    },
    'minimalprojectionarea': {
        'type': 'num',
        'name': 'Minimal Projection Area',
        'calculatorSoftware': calculatorSoftware,
        'calculatorSoftwareVersion': '15_6',
        'minimum': 0,
    },
    'minimalprojectionradius': {
        'type': 'num',
        'name': 'Minimal Projection Radius',
        'calculatorSoftware': calculatorSoftware,
        'calculatorSoftwareVersion': '15_6',
        'minimum': 0,
    },
    'minimalprojectionsize': {
        'type': 'num',
        'name': 'Minimal Projection Size',
        'calculatorSoftware': calculatorSoftware,
        'calculatorSoftwareVersion': '15_6',
        'minimum': 0,
    }
}

# The descriptors listed below also have the option to be calculated
# at a variety of pH values, which yields a distribution.
# This will be done, and they will be added as different descriptors
_pHDependentDescriptors = {
    'avgpol': {
        'type': 'num',
        'name': 'Average Molecular Polarizability',
        'calculatorSoftware': calculatorSoftware,
        'calculatorSoftwareVersion': '15_6',
        'maximum': None,
        'minimum': None,
    },
    'molpol': {
        'type': 'num',
        'name': 'Specific Molecular Polarizability',
        'calculatorSoftware': calculatorSoftware,
        'calculatorSoftwareVersion': '15_6',
        'maximum': None,
        'minimum': None,
    },
    'vanderwaals': {
        'type': 'num',
        'name': 'Van der Waals Surface Area',
        'calculatorSoftware': calculatorSoftware,
        'calculatorSoftwareVersion': '15_6',
        'maximum': None,
        'minimum': 0,
    },
    'asa': {
        'type': 'num',
        'name': 'Water Acessible Surface Area',
        'calculatorSoftware': calculatorSoftware,
        'calculatorSoftwareVersion': '15_6',
        'maximum': None,
        'minimum': 0,
    },
    'asa+': {
        'type': 'num',
        'name': 'Partial Positive Charged water accessible surface area',
        'calculatorSoftware': calculatorSoftware,
        'calculatorSoftwareVersion': '15_6',
        'maximum': None,
        'minimum': 0,
    },
    'asa-': {
        'type': 'num',
        'name': 'Partial negative Charged water accessible surface area',
        'calculatorSoftware': calculatorSoftware,
        'calculatorSoftwareVersion': '15_6',
        'maximum': None,
        'minimum': 0,
    },
    'asa_hydrophobic': {
        'type': 'num',
        'name': 'Hydrophobic water accessible surface area',
        'calculatorSoftware': calculatorSoftware,
        'calculatorSoftwareVersion': '15_6',
        'maximum': None,
        'minimum': 0,
    },
    'asa_polar': {
        'type': 'num',
        'name': 'Polar water accessible surface area',
        'calculatorSoftware': calculatorSoftware,
        'calculatorSoftwareVersion': '15_6',
        'maximum': None,
        'minimum': 0,
    },
    'hbda_acc': {
        'type': 'num',
        'name': 'Hydrogen bond acceptor count',
        'calculatorSoftware': calculatorSoftware,
        'calculatorSoftwareVersion': '15_6',
        'maximum': None,
        'minimum': 0,
    },
    'hbda_don': {
        'type': 'num',
        'name': 'Hydrogen bond donor count',
        'calculatorSoftware': calculatorSoftware,
        'calculatorSoftwareVersion': '15_6',
        'maximum': None,
        'minimum': 0,
    },
    'polar_surface_area': {
        'type': 'num',
        'name': 'Polar Surface Area',
        'calculatorSoftware': calculatorSoftware,
        # this is not a typo. This descriptor was introduced later
        'calculatorSoftwareVersion': '16_5',
        'maximum': None,
        'minimum': 0,
    }
}

cxcalcCommands = OrderedDict()

for key in _descriptorDict.keys():
    cxcalcCommands[key] = key


def setup_pHdependentDescriptors(_descriptorDict):
    """Set up for calculation of pH dependent descriptors."""
    pH_vals = DRP.models.NumRxnDescriptorValue.objects.filter(descriptor__heading='reaction_pH', reaction__performedreaction__valid=True).exclude(
        value=None).order_by('value').values_list('value', flat=True).distinct()
    for descriptor, d in _pHDependentDescriptors.items():
        for pH in pH_vals:
            pH_string = str(pH).replace('.', '_')  # R compatibility
            d_copy = d.copy()
            d_copy['name'] += ' at pH {}'.format(pH_string)
            _descriptorDict[descriptor + '_pH{}'.format(pH_string)] = d_copy

    descriptorDict = setup(_descriptorDict)

    _cxcalcpHCommandStems = {
        'avgpol_pH{}': 'avgpol -H {}',
        'molpol_pH{}': 'molpol -H {}',
        'vanderwaals_pH{}': 'vdwsa -H {}',
        'asa_pH{}': 'molecularsurfacearea -t ASA -H {}',
        'asa+_pH{}': 'molecularsurfacearea -t ASA+ -H {}',
        'asa-_pH{}': 'molecularsurfacearea -t ASA- -H {}',
        'asa_hydrophobic_pH{}': 'molecularsurfacearea -t ASA_H -H {}',
        'asa_polar_pH{}': 'molecularsurfacearea -t ASA_P -H {}',
        'hbda_acc_pH{}': 'acceptorcount -H {}',
        'hbda_don_pH{}': 'donorcount -H {}',
        'polar_surface_area_pH{}': 'polarsurfacearea -H {}',
    }

    for key, command in _cxcalcpHCommandStems.items():
        for pH in pH_vals:
            pH_string = str(pH).replace('.', '_')  # R compatibility
            cxcalcCommands[key.format(pH_string)] = command.format(pH)

    if len(cxcalcCommands) != len(_descriptorDict):
        raise RuntimeError(
            "Need the same number of cxcalc commands as descriptors being calculated")

    return descriptorDict


def delete_descriptors(compound_set, descriptorDict, cxcalcCommands):
    """Bulk deletion of descriptors."""
    DRP.models.NumMolDescriptorValue.objects.filter(descriptor__in=[descriptorDict[ck] for ck in cxcalcCommands.keys() if _descriptorDict[ck]['type'] == 'num'],
                                                    compound__in=compound_set).delete(recalculate_reactions=False)
    DRP.models.OrdMolDescriptorValue.objects.filter(descriptor__in=[descriptorDict[ck] for ck in cxcalcCommands.keys() if _descriptorDict[ck]['type'] == 'ord'],
                                                    compound__in=compound_set).delete(recalculate_reactions=False)
    DRP.models.BoolMolDescriptorValue.objects.filter(descriptor__in=[descriptorDict[ck] for ck in cxcalcCommands.keys() if _descriptorDict[ck]['type'] == 'bool'],
                                                     compound__in=compound_set).delete(recalculate_reactions=False)


def calculate_many(compound_set, verbose=False, whitelist=None):
    """Bulk calculation of descriptors."""
    if verbose:
        print "Creating descriptor dictionary"
    descriptorDict = setup_pHdependentDescriptors(_descriptorDict)
    if whitelist is not None:
        filtered_cxcalcCommands = {k: cxcalcCommands[
            k] for k in cxcalcCommands.keys() if k in whitelist}
    else:
        filtered_cxcalcCommands = cxcalcCommands
    if verbose:
        print "Deleting old descriptor values."
    delete_descriptors(compound_set, descriptorDict, cxcalcCommands)

    num_to_create = []
    ord_to_create = []
    for i, compound in enumerate(compound_set):
        if verbose:
            print "{}; Compound {} ({}/{})".format(compound, compound.pk, i + 1, len(compound_set))
        num_to_create, ord_to_create = _calculate(
            compound, descriptorDict, filtered_cxcalcCommands, verbose=verbose, num_to_create=num_to_create, ord_to_create=ord_to_create)
        if len(num_to_create) > create_threshold:
            if verbose:
                print 'Creating {} numeric values'.format(len(num_to_create))
            DRP.models.NumMolDescriptorValue.objects.bulk_create(num_to_create)
            num_to_create = []
        if len(ord_to_create) > create_threshold:
            if verbose:
                print 'Creating {} ordinal values'.format(len(Ord_to_create))
            DRP.models.OrdMolDescriptorValue.objects.bulk_create(ord_to_create)
            ord_to_create = []

    if verbose:
        print 'Creating {} numeric values'.format(len(num_to_create))
    DRP.models.NumMolDescriptorValue.objects.bulk_create(num_to_create)
    if verbose:
        print 'Creating {} ordinal values'.format(len(ord_to_create))
    DRP.models.OrdMolDescriptorValue.objects.bulk_create(ord_to_create)


def calculate(compound, verbose=False, whitelist=None):
    """Calculate descriptor values."""
    if verbose:
        print "Creating descriptor dictionary"
    descriptorDict = setup_pHdependentDescriptors(_descriptorDict)
    if whitelist is not None:
        filtered_cxcalcCommands = {k: cxcalcCommands[
            k] for k in cxcalcCommands.keys() if k in whitelist}
    else:
        filtered_cxcalcCommands = cxcalcCommands
    if verbose:
        print "Deleting old descriptor values"
    delete_descriptors([compound], descriptorDict, cxcalcCommands)
    if verbose:
        print "Creating new descriptor values."
    num_to_create, ord_to_create = _calculate(
        compound, descriptorDict, filtered_cxcalcCommands, verbose=verbose)

    if verbose:
        print "Creating {} numerical and {} ordinal".format(len(num_to_create), len(ord_to_create))
    DRP.models.NumMolDescriptorValue.objects.bulk_create(num_to_create)
    DRP.models.OrdMolDescriptorValue.objects.bulk_create(ord_to_create)


def _calculate(compound, descriptorDict, cxcalcCommands, verbose=False, num_to_create=[], ord_to_create=[]):
    notFound = True
    if notFound and (compound.smiles is not None and compound.smiles != ''):
        lecProc = Popen([settings.CHEMAXON_DIR[CHEMAXON_VERSION] + 'cxcalc', compound.smiles,
                         'leconformer'], stdout=PIPE, stderr=PIPE, close_fds=True)  # lec = lowest energy conformer
        lecProc.wait()
        if lecProc.returncode == 0:
            lec, lecErr = lecProc.communicate()
            notFound = False
    if notFound and (compound.INCHI is not None and compound.INCHI != ''):
        lecProc = Popen([settings.CHEMAXON_DIR[CHEMAXON_VERSION] + 'cxcalc', compound.INCHI,
                         'leconformer'], stdout=PIPE, stderr=PIPE, close_fds=True)  # lec = lowest energy conformer
        lecProc.wait()
        if lecProc.returncode == 0:
            lec, lecErr = lecProc.communicate()
            notFound = False
    if not notFound:
        # -N ih means leave off the header row and id column
        calcProc = Popen([settings.CHEMAXON_DIR[CHEMAXON_VERSION] + 'cxcalc', '-N', 'ih', lec] + [x for x in chain(
            *(command.split(' ') for command in cxcalcCommands.values()))], stdout=PIPE, stderr=PIPE, close_fds=True)
        calcProc.wait()
        if calcProc.returncode == 0:
            res, resErr = calcProc.communicate()
            if not resErr:
                resLines = res.split('\n')
                if len(resLines) == 2:  # last line is blank
                    resList = resLines[0].split('\t')
                    commandKeys = cxcalcCommands.keys()

                    if len(resList) == len(commandKeys):
                        for i in range(len(resList)):
                            if _descriptorDict[commandKeys[i]]['type'] == 'num':
                                n = DRP.models.NumMolDescriptorValue(descriptor=descriptorDict[commandKeys[
                                                                     i]], compound=compound, value=float(resList[i]))
                                try:
                                    n.full_clean()
                                except ValidationError as e:
                                    warnings.warn('Value {} for compound {} and descriptor {} failed validation. Value set to None. Validation error message: {}'.format(
                                        n.value, n.compound, n.descriptor, e))
                                    n.value = None
                                num_to_create.append(n)
                            elif _descriptorDict[commandKeys[i]]['type'] == 'ord':
                                o = DRP.models.OrdMolDescriptorValue(descriptor=descriptorDict[commandKeys[
                                                                     i]], compound=compound, value=int(resList[i]))
                                try:
                                    o.full_clean()
                                except ValidationError as e:
                                    warnings.warn('Value {} for compound {} and descriptor {} failed validation. Value set to None. Validation error message: {}'.format(
                                        o.value, o.compound, o.descriptor, e))
                                    o.value = None
                                ord_to_create.append(o)
                            else:
                                return ValueError('Descriptor has unrecognized type {}'.format(_descriptorDict[commandKeys[i]]['type']))
                            # elif _descriptorDict[commandKeys[i]]['type'] == 'bool':
                                # TODO Not sure whether cxcalc even returns any boolean values, but if it does I don't know how it notates them and they should be coerced correctly
                                # commenting out this bit since it should be double checked before anyone uses it
                                # bool_to_create.append(DRP.models.BoolMolDescriptorValue(descriptor=descriptorDict[commandKeys[i]], compound=compound, bool(int(value=resList[i])))
                            # NOTE: No categorical descriptors are included yet, and since they are more complicated to code I've left it for the moment.
                            # NOTE: Calculation failure values are not included in the documentation, so I've assumed that it doesn't happen, since we have no way of identifying
                            # for it other than for the database to push it out
                            # as a part of validation procedures.

                    else:
                        raise RuntimeError("Number of cxcalc commands ({}) does not match number of results ({})".format(
                            len(commandKeys), len(resList)))
            else:
                warnings.warn("cxcalc returned error: {}".format(resErr))
        else:
            warnings.warn("cxcalc exited with nonzero return code {}".format(
                calcProc.returncode))
    else:
        warnings.warn("Compound not found")

    return num_to_create, ord_to_create
