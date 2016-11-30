"""A descriptor calculator to wrap and use calls to ChemAxon.

Requires cxcalc (part of JChem) to be installed and licensed.

"""
import DRP
from .utils import setup
from django.conf import settings
from django.core.exceptions import ValidationError
from collections import OrderedDict
from subprocess import Popen, PIPE
from itertools import chain
import logging
logger = logging.getLogger("DRP")

# The version of ChemAxon currently in use.
# Should be matched by a entry in the dictionary in the settings file
# This version is reflected in, but not necessarily the same as, the
# descriptor versions
MIN_CHEMAXON_VERSION = '16.5'
MAX_CHEMAXON_VERSION = '16.10'

MIN_CHEMAXON_VERSION = tuple(int(element) for element in MIN_CHEMAXON_VERSION.split('.'))
MAX_CHEMAXON_VERSION = tuple(int(element) for element in MAX_CHEMAXON_VERSION.split('.'))

CHEMAXON_VERSION = None
for key in settings.CHEMAXON_DIR.keys():
    tuplekey = tuple(int(element) for element in key.split('.'))
    if tuplekey >= MIN_CHEMAXON_VERSION and tuplekey <= MAX_CHEMAXON_VERSION:
        if CHEMAXON_VERSION == None or tuple(int(element) for element in CHEMAXON_VERSION.split('.')) < tuplekey:
             CHEMAXON_VERSION = key

if CHEMAXON_VERSION == None:
    raise ValueError('There must be a chemaxon install with versions between {} and {}'.format(MIN_CHEMAXON_VERSION, MAX_CHEMAXON_VERSION))


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

_descriptorDict['vdw_area_N_ratio'] = {
    'type': 'num',
    'name': 'vdw area divided by nitrogen count',
    'calculatorSoftware': 'drp_chemaxon',
    'calculatorSoftwareVersion':'16_5',
    'maximum':None,
    'minimum':0
}

def setup_pHdependentDescriptors(_descriptorDict):
    """Set up for calculation of pH dependent descriptors."""
    pH_vals = DRP.models.NumRxnDescriptorValue.objects.filter(descriptor__heading='reaction_pH', reaction__performedreaction__valid=True).exclude(
        value=None).order_by('value').values_list('value', flat=True).distinct()
    for descriptor, d in _pHDependentDescriptors.items():
        d_copy = d.copy()
        d_copy['name'] += ' for nominal structure'
        _descriptorDict[descriptor + '_nominal'] = d_copy

        for pH in pH_vals:
            pH_string = str(pH).replace('.', '_')  # R compatibility
            d_copy = d.copy()
            d_copy['name'] += ' at pH {}'.format(pH_string)
            _descriptorDict[descriptor + '_pH{}'.format(pH_string)] = d_copy

    descriptorDict = setup(_descriptorDict)

    _cxcalcpHCommandStems = {
        'avgpol': 'avgpol',
        'molpol': 'molpol',
        'vanderwaals': 'vdwsa',
        'asa': 'molecularsurfacearea -t ASA',
        'asa+': 'molecularsurfacearea -t ASA+',
        'asa-': 'molecularsurfacearea -t ASA-',
        'asa_hydrophobic': 'molecularsurfacearea -t ASA_H',
        'asa_polar': 'molecularsurfacearea -t ASA_P',
        'hbda_acc': 'acceptorcount',
        'hbda_don': 'donorcount',
        'polar_surface_area': 'polarsurfacearea',
    }

    for key, command in _cxcalcpHCommandStems.items():
        # nominal structure
        cxcalcCommands["{}_nominal".format(key)] = command

        for pH in pH_vals:
            pH_string = str(pH).replace('.', '_')  # R compatibility
            cxcalcCommands["{}_pH{}".format(key, pH_string)] = "{} -H {}".format(command, pH)

    return descriptorDict


def delete_descriptors(compound_set, descriptorDict):
    """Bulk deletion of descriptors."""
    DRP.models.NumMolDescriptorValue.objects.filter(descriptor__in=[descriptorDict[ck] for ck in descriptorDict if _descriptorDict[ck]['type'] == 'num'],
                                                    compound__in=compound_set).delete()
    DRP.models.OrdMolDescriptorValue.objects.filter(descriptor__in=[descriptorDict[ck] for ck in descriptorDict if  _descriptorDict[ck]['type'] == 'ord'],
                                                    compound__in=compound_set).delete()
    DRP.models.BoolMolDescriptorValue.objects.filter(descriptor__in=[descriptorDict[ck] for ck in descriptorDict if _descriptorDict[ck]['type'] == 'bool'],
                                                     compound__in=compound_set).delete()


def calculate_many(compound_set, verbose=False, whitelist=None):
    """Bulk calculation of descriptors."""
    if verbose:
        logger.info("Creating descriptor dictionary")
    descriptorDict = setup_pHdependentDescriptors(_descriptorDict)
    if whitelist is None:
        descriptorDict = {k:v for k,v in descriptorDict.items()}
    else:
        descriptorDict = {k:v for k,v in descriptorDict.items() if k in whitelist}
    if verbose:
        logger.info("Deleting old descriptor values.")
    delete_descriptors(compound_set, descriptorDict)

    num_to_create = []
    ord_to_create = []
    filtered_cxcalcCommands = {k:v for k,v in cxcalcCommands.items() if k in descriptorDict.keys()}
    for i, compound in enumerate(compound_set):
        if verbose:
            logger.info("{}; Compound {} ({}/{})".format(compound, compound.pk, i + 1, len(compound_set)))
        num_to_create, ord_to_create = _calculate(
            compound, descriptorDict, filtered_cxcalcCommands, verbose=verbose, num_to_create=num_to_create, ord_to_create=ord_to_create)
        if len(num_to_create) > create_threshold:
            if verbose:
                logger.info('Creating {} numeric values'.format(len(num_to_create)))
            DRP.models.NumMolDescriptorValue.objects.bulk_create(num_to_create)
            num_to_create = []
        if len(ord_to_create) > create_threshold:
            if verbose:
                logger.info('Creating {} ordinal values'.format(len(Ord_to_create)))
            DRP.models.OrdMolDescriptorValue.objects.bulk_create(ord_to_create)
            ord_to_create = []

    if verbose:
        logger.info('Creating {} numeric values'.format(len(num_to_create)))
    DRP.models.NumMolDescriptorValue.objects.bulk_create(num_to_create)
    if verbose:
        logger.info('Creating {} ordinal values'.format(len(ord_to_create)))
    DRP.models.OrdMolDescriptorValue.objects.bulk_create(ord_to_create)


def calculate(compound, verbose=False, whitelist=None):
    """Calculate descriptor values."""
    if verbose:
        logger.info("Creating descriptor dictionary")
    descriptorDict = setup_pHdependentDescriptors(_descriptorDict)
    if whitelist is None:
        descriptorDict = {k:v for k,v in descriptorDict.items()}
    else:
        descriptorDict = {k:v for k,v in descriptorDict.items() if k in whitelist}
    if verbose:
        logger.info("Deleting old descriptor values.")
    delete_descriptors([compound], descriptorDict, cxcalcCommands)
    filtered_cxcalcCommands = {k:v for k,v in cxcalcCommands.items() if k in descriptorDict.keys()}
    if verbose:
        logger.info("Creating new descriptor values.")
    num_to_create, ord_to_create = _calculate(
        compound, descriptorDict, filtered_cxcalcCommands, verbose=verbose)

    if verbose:
        logger.info("Creating {} numerical and {} ordinal".format(len(num_to_create), len(ord_to_create)))
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
    if not notFound or lec != '':
        # -N ih means leave off the header row and id column
        calcProc = Popen([settings.CHEMAXON_DIR[CHEMAXON_VERSION] + 'cxcalc', '-N', 'ih', lec] + [x for x in chain(
            *(command.split(' ') for command in cxcalcCommands.values()))], stdout=PIPE, stderr=PIPE, close_fds=True)
        calcProc.wait()
        if calcProc.returncode == 0:
            res, resErr = calcProc.communicate()
            if not resErr:
                resLines = res.decode('UTF-8').split('\n')
                if len(resLines) == 2:  # last line is blank
                    resList = resLines[0].split('\t')
                    commandKeys = tuple(cxcalcCommands.keys())

                    if len(resList) == len(commandKeys):
                        for i in range(len(resList)):
                            if _descriptorDict[commandKeys[i]]['type'] == 'num':
                                n = DRP.models.NumMolDescriptorValue(descriptor=descriptorDict[commandKeys[
                                                                     i]], compound=compound, value=float(resList[i]))
                                if commandKeys[i] == 'vanderwaals' and 'N' in compound.elements.keys(): # I hate this special case, but this might not stick around so I'm leaving it for now
                                    n2 = DRP.models.NumMolDescriptorValue(descriptor=descriptorDict['vdw_area_N_ratio'], compound=compound, 
                                                                     value= float(resList[i])/compound.elements['N']['stoichiometry'])  
                                else:
                                    n2=None
                                try:
                                    n.full_clean()
                                    if n2 is not None:
                                        n2.full_clean()
                                except ValidationError as e:
                                    logger.warning('Value {} for compound {} and descriptor {} failed validation. Value set to None. Validation error message: {}'.format(
                                        n.value, n.compound, n.descriptor, e))
                                    n.value = None
                                num_to_create.append(n)
                                if n2 is not None:
                                    num_to_create.append(n2)
                            elif _descriptorDict[commandKeys[i]]['type'] == 'ord':
                                o = DRP.models.OrdMolDescriptorValue(descriptor=descriptorDict[commandKeys[
                                                                     i]], compound=compound, value=int(resList[i]))
                                try:
                                    o.full_clean()
                                except ValidationError as e:
                                    logger.warning('Value {} for compound {} and descriptor {} failed validation. Value set to None. Validation error message: {}'.format(
                                        o.value, o.compound, o.descriptor, e))
                                    o.value = None
                                ord_to_create.append(o)
                            else:
                                return ValueError('Descriptor has unrecognized type {}'.format(_descriptorDict[commandKeys[i]]['type']))
                            # elif _descriptorDict[commandKeys[i]]['type'] == 'bool':
                                # Not sure whether cxcalc even returns any boolean values, but if it does I don't know how it notates them and they should be coerced correctly
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
                logger.warning("cxcalc returned error: {}".format(resErr))
        else:
            logger.warning("cxcalc exited with nonzero return code {}".format(
                calcProc.returncode))
    else:
        logger.warning("Compound not found")

    return num_to_create, ord_to_create
