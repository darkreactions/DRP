"""A descriptor calculator to wrap and use calls to ChemAxon.

Requires cxcalc (part of JChem) to be installed and licensed.

"""
import DRP
from utils import setup
from django.conf import settings
from collections import OrderedDict
from subprocess import Popen, PIPE
from itertools import chain

_descriptorDict = {
    'refractivity': {
        'type': 'num',
        'name': 'Refractivity',
        'calculatorSoftware': 'ChemAxon/cxcalc',
        'calculatorSoftwareVersion': '15.6',
    },
    'maximalprojectionarea': {
        'type': 'num',
        'name': 'Maximal Projection Area',
        'calculatorSoftware': 'ChemAxon/cxcalc',
        'calculatorSoftwareVersion': '15.6'
    },
    'maximalprojectionradius': {
        'type': 'num',
        'name': 'Maximal Projection Radius',
        'calculatorSoftware': 'ChemAxon/cxcalc',
        'calculatorSoftwareVersion': '15.6'
    },
    'maximalprojectionsize': {
        'type': 'num',
        'name': 'Maximal Projection Size',
        'calculatorSoftware': 'ChemAxon/cxcalc',
        'calculatorSoftwareVersion': '15.6'
    },
    'minimalprojectionarea': {
        'type': 'num',
        'name': 'Minimal Projection Area',
        'calculatorSoftware': 'ChemAxon/cxcalc',
        'calculatorSoftwareVersion': '15.6'
    },
    'minimalprojectionradius': {
        'type': 'num',
        'name': 'Minimal Projection Radius',
        'calculatorSoftware': 'ChemAxon/cxcalc',
        'calculatorSoftwareVersion': '15.6'
    },
    'minimalprojectionsize': {
        'type': 'num',
        'name': 'Minimal Projection Size',
        'calculatorSoftware': 'ChemAxon/cxcalc',
        'calculatorSoftwareVersion': '15.6'
    }
}

#The descriptors listed below also have the option to be calculated
#at a variety of pH values, which yields a distribution.
#This will be done, and they will be added as different descriptors
_pHDependantDescriptors = {
    'avgpol': {
        'type': 'num',
        'name': 'Average Molecular Polarizability',
        'calculatorSoftware': 'ChemAxon/cxcalc',
        'calculatorSoftwareVersion': '15.6'
    },
    'molpol': {
        'type': 'num',
        'name': 'Specific Molecular Polarizability',
        'calculatorSoftware': 'ChemAxon/cxcalc',
        'calculatorSoftwareVersion': '15.6'
    },
    'vanderwaals': {
        'type': 'num',
        'name': 'Van der Waals Surface Area',
        'calculatorSoftware': 'ChemAxon/cxcalc',
        'calculatorSoftwareVersion': '15.6'
    },
    'asa': {
        'type': 'num',
        'name': 'Water Acessible Surface Area',
        'calculatorSoftware': 'ChemAxon/cxcalc',
        'calculatorSoftwareVersion': '15.6'
    },
    'asa+': {
        'type': 'num',
        'name': 'Partial Positive Charged water accessible surface area',
        'calculatorSoftware': 'ChemAxon/cxcalc',
        'calculatorSoftwareVersion': '15.6'
    },
    'asa-': {
        'type': 'num',
        'name': 'Partial negative Charged water accessible surface area',
        'calculatorSoftware': 'ChemAxon/cxcalc',
        'calculatorSoftwareVersion': '15.6'
    },
    'asa_hydrophobic': {
        'type': 'num',
        'name': 'Hydrophobic water accessible surface area',
        'calculatorSoftware': 'ChemAxon/cxcalc',
        'calculatorSoftwareVersion': '15.6'
    },
    'asa_polar': {
        'type': 'num',
        'name': 'Polar water accessible surface area',
        'calculatorSoftware': 'ChemAxon/cxcalc',
        'calculatorSoftwareVersion': '15.6'
    },
    'hbda_acc': {
        'type': 'num',
        'name': 'Hydrogen bond acceptor count',
        'calculatorSoftware': 'ChemAxon/cxcalc',
        'calculatorSoftwareVersion': '15.6'
    },
    'hbda_don': {
        'type': 'num',
        'name': 'Hydrogen bond donor count',
        'calculatorSoftware': 'ChemAxon/cxcalc',
        'calculatorSoftwareVersion': '15.6'
    }
}


for descriptor, d in _pHDependantDescriptors.items():
    for i in range(1, 15):
        d['name'] += ' at pH {}'.format(i)
        _descriptorDict[descriptor + '_pH{}'.format(i)] = d

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
    'hbda_don_pH{}': 'donorcount -H {}'
}

cxcalcCommands = OrderedDict()

for key, command in _cxcalcpHCommandStems.items():
    for i in range(1, 15):
        cxcalcCommands[key.format(i)] = command.format(i)


def calculate(compound):
    notFound = True
    if notFound and (compound.INCHI is not None and compound.INCHI is not ''):
        lecProc = Popen([settings.CHEMAXON_DIR['15.6'] + 'cxcalc', compound.INCHI, 'leconformer'], stdout=PIPE, close_fds=True) # lec = lowest energy conformer
        lecProc.wait()
        if lecProc.returncode == 0:
            lec, lecErr = lecProc.communicate()
            notFound = False
    if notFound and (compound.smiles is not None and compound.smiles is not ''):
        lecProc = Popen([settings.CHEMAXON_DIR['15.6'] + 'cxcalc', compound.smiles, 'leconformer'], stdout=PIPE, close_fds=True) # lec = lowest energy conformer
        lecProc.wait()
        if lecProc.returncode == 0:
            lec, lecErr = lecProc.communicate()
            notFound = False  
        notFound = False
    if not notFound:
        calcProc = Popen([settings.CHEMAXON_DIR['15.6'] + 'cxcalc', lec] + [x for x in chain(*(command.split(' ') for command in cxcalcCommands.values()))], stdout=PIPE, close_fds=True) 
        calcProc.wait()
        if calcProc.returncode == 0:
            res, resErr = calcProc.communicate()
            if resErr == '':
                resLines = res.split('\n')
                if len(resLines) == 3:
                    resList = resLines[1].split('\t')
                    commandKeys = cxcalcCommands.keys()
                    for i in len(resList):
                        if _descriptorDict[commandKeys[i]]['type'] == 'num':
                            DRP.models.NumMolDescriptorValue.get_or_create(descriptor=descriptorDict[commandKeys[i]], compound=compound, value=resList[i])
                        elif _descriptorDict[commandKeys[i]]['type'] == 'ord':
                            DRP.models.OrdMolDescriptorValue.get_or_create(descriptor=descriptorDict[commandKeys[i]], compound=compound, value=resList[i])
                        elif _descriptorDict[commandKeys[i]]['type'] == 'bool':
                            DRP.models.BoolMolDescriptorValue.get_or_create(descriptor=descriptorDict[commandKeys[i]], compound=compound, value=resList[i])
                        # NOTE: No categorical descriptors are included yet, and since they are more complicated to code I've left it for the moment.
                        # NOTE: Calculation failure values are not included in the documentation, so I've assumed that it doesn't happen, since we have no way of identifying
                        # for it other than for the database to push it out as a part of validation procedures.
