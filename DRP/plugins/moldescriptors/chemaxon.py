"""A descriptor calculator to wrap and use calls to ChemAxon.

Requires cxcalc (part of JChem) to be installed and licensed.

"""
import DRP
from utils import setup
from django.conf import settings
from collections import OrderedDict
from subprocess import Popen, PIPE
from itertools import chain
import warnings

_descriptorDict = {
    'refractivity': {
        'type': 'num',
        'name': 'Refractivity',
        'calculatorSoftware': 'ChemAxon_cxcalc',
        'calculatorSoftwareVersion': '15.6',
    },
    'maximalprojectionarea': {
        'type': 'num',
        'name': 'Maximal Projection Area',
        'calculatorSoftware': 'ChemAxon_cxcalc',
        'calculatorSoftwareVersion': '15.6'
    },
    'maximalprojectionradius': {
        'type': 'num',
        'name': 'Maximal Projection Radius',
        'calculatorSoftware': 'ChemAxon_cxcalc',
        'calculatorSoftwareVersion': '15.6'
    },
    'maximalprojectionsize': {
        'type': 'num',
        'name': 'Maximal Projection Size',
        'calculatorSoftware': 'ChemAxon_cxcalc',
        'calculatorSoftwareVersion': '15.6'
    },
    'minimalprojectionarea': {
        'type': 'num',
        'name': 'Minimal Projection Area',
        'calculatorSoftware': 'ChemAxon_cxcalc',
        'calculatorSoftwareVersion': '15.6'
    },
    'minimalprojectionradius': {
        'type': 'num',
        'name': 'Minimal Projection Radius',
        'calculatorSoftware': 'ChemAxon_cxcalc',
        'calculatorSoftwareVersion': '15.6'
    },
    'minimalprojectionsize': {
        'type': 'num',
        'name': 'Minimal Projection Size',
        'calculatorSoftware': 'ChemAxon_cxcalc',
        'calculatorSoftwareVersion': '15.6'
    }
}

#The descriptors listed below also have the option to be calculated
#at a variety of pH values, which yields a distribution.
#This will be done, and they will be added as different descriptors
_pHDependentDescriptors = {
    'avgpol': {
        'type': 'num',
        'name': 'Average Molecular Polarizability',
        'calculatorSoftware': 'ChemAxon_cxcalc',
        'calculatorSoftwareVersion': '15.6'
    },
    'molpol': {
        'type': 'num',
        'name': 'Specific Molecular Polarizability',
        'calculatorSoftware': 'ChemAxon_cxcalc',
        'calculatorSoftwareVersion': '15.6'
    },
    'vanderwaals': {
        'type': 'num',
        'name': 'Van der Waals Surface Area',
        'calculatorSoftware': 'ChemAxon_cxcalc',
        'calculatorSoftwareVersion': '15.6'
    },
    'asa': {
        'type': 'num',
        'name': 'Water Acessible Surface Area',
        'calculatorSoftware': 'ChemAxon_cxcalc',
        'calculatorSoftwareVersion': '15.6'
    },
    'asa+': {
        'type': 'num',
        'name': 'Partial Positive Charged water accessible surface area',
        'calculatorSoftware': 'ChemAxon_cxcalc',
        'calculatorSoftwareVersion': '15.6'
    },
    'asa-': {
        'type': 'num',
        'name': 'Partial negative Charged water accessible surface area',
        'calculatorSoftware': 'ChemAxon_cxcalc',
        'calculatorSoftwareVersion': '15.6'
    },
    'asa_hydrophobic': {
        'type': 'num',
        'name': 'Hydrophobic water accessible surface area',
        'calculatorSoftware': 'ChemAxon_cxcalc',
        'calculatorSoftwareVersion': '15.6'
    },
    'asa_polar': {
        'type': 'num',
        'name': 'Polar water accessible surface area',
        'calculatorSoftware': 'ChemAxon_cxcalc',
        'calculatorSoftwareVersion': '15.6'
    },
    'hbda_acc': {
        'type': 'num',
        'name': 'Hydrogen bond acceptor count',
        'calculatorSoftware': 'ChemAxon_cxcalc',
        'calculatorSoftwareVersion': '15.6'
    },
    'hbda_don': {
        'type': 'num',
        'name': 'Hydrogen bond donor count',
        'calculatorSoftware': 'ChemAxon_cxcalc',
        'calculatorSoftwareVersion': '15.6'
    }
}

cxcalcCommands = OrderedDict()

for key in _descriptorDict.keys():
    cxcalcCommands[key] = key

for descriptor, d in _pHDependentDescriptors.items():
    for i in range(1, 15):
        d_copy = d.copy()
        d_copy['name'] += ' at pH {}'.format(i)
        _descriptorDict[descriptor + '_pH{}'.format(i)] = d_copy

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


for key, command in _cxcalcpHCommandStems.items():
    for i in range(1, 15):
        cxcalcCommands[key.format(i)] = command.format(i)

if len(cxcalcCommands) != len(_descriptorDict):
    raise RuntimeError("Need the same number of cxcalc commands as descriptors being calculated")

def delete_descriptors(compound_set):
    DRP.models.NumMolDescriptorValue.objects.filter(descriptor__in=[descriptorDict[ck] for ck in cxcalcCommands.keys() if _descriptorDict[ck]['type']=='num'],
                                                    compound__in=compound_set).delete(recalculate_reactions=False)
    DRP.models.OrdMolDescriptorValue.objects.filter(descriptor__in=[descriptorDict[ck] for ck in cxcalcCommands.keys() if _descriptorDict[ck]['type']=='ord'],
                                                    compound__in=compound_set).delete(recalculate_reactions=False)
    DRP.models.BoolMolDescriptorValue.objects.filter(descriptor__in=[descriptorDict[ck] for ck in cxcalcCommands.keys() if _descriptorDict[ck]['type']=='bool'],
                                                     compound__in=compound_set).delete(recalculate_reactions=False)

def calculate_many(compound_set, verbose=False):
    if verbose:
        print "Deleting old descriptors"
    delete_descriptors(compound_set)
    # filtering first and then doing a big cxcalc command is actually slower
    #if verbose:
        #print "Filtering out compounds without valid SMILES or INCHI"
    #filtered_compounds = _filter_compounds(compound_set, verbose=verbose)
    #_calculate_with_lec(filtered_compounds, verbose=verbose)
    for i, compound in enumerate(compound_set):
        if verbose:
            print "{}; Compound {} ({}/{})".format(compound, compound.pk, i+1, len(compound_set))
        _calculate(compound, verbose=verbose)

def calculate(compound):
    delete_descriptors([compound])
    _calculate(compound)

def _filter_compounds(compound_set, verbose=False):
    filtered_compounds = OrderedDict()
    for i, compound in enumerate(compound_set):
        notFound = True
        if notFound and (compound.smiles is not None and compound.smiles != ''):
            lecProc = Popen([settings.CHEMAXON_DIR['15.6'] + 'cxcalc', compound.smiles, 'leconformer'], stdout=PIPE, stderr=PIPE, close_fds=True) # lec = lowest energy conformer
            lecProc.wait()
            if lecProc.returncode == 0:
                lec, lecErr = lecProc.communicate()
                notFound = False  
        if notFound and (compound.INCHI is not None and compound.INCHI != ''):
            lecProc = Popen([settings.CHEMAXON_DIR['15.6'] + 'cxcalc', compound.INCHI, 'leconformer'], stdout=PIPE, stderr=PIPE, close_fds=True) # lec = lowest energy conformer
            lecProc.wait()
            if lecProc.returncode == 0:
                lec, lecErr = lecProc.communicate()
                notFound = False
        if not notFound:
            filtered_compounds[lec] = compound
        else:
            warnings.warn("Compound not found")
        if verbose:
            print "Done with {} ({}/{})".format(compound, i+1, len(compound_set))
    return filtered_compounds


def _calculate_with_lec(compound_dict, verbose=False):
    """
    Does cxcalc calculations given an ordered dictionary lec representation : compound
    """
    if verbose:
        print "Calculating descriptors for {} compounds with valid SMILES/INCHI".format(len(compound_dict))
    #-N ih says not to display id or header
    calcProc = Popen([settings.CHEMAXON_DIR['15.6'] + 'cxcalc', '-N', 'ih'] + compound_dict.keys() +
                     [x for x in chain(*(command.split(' ') for command in cxcalcCommands.values()))], stdout=PIPE, stderr=PIPE, close_fds=True)
    calcProc.wait()
    if calcProc.returncode == 0:
        res, resErr = calcProc.communicate()
        if not resErr:
            if verbose:
                print "Putting descriptor values in database"
            resLines = res.split('\n')[:-1] #the last line is blank
            if len(resLines) == len(compound_dict): 
                commandKeys = cxcalcCommands.keys()
                for i, resLine in enumerate(resLines): 
                    resList = resLine.split('\t')
                    compound = compound_dict.values()[i]
    
                    num_to_create = []
                    ord_to_create = []
                    bool_to_create = []
                    for i in range(len(resList)):
                        if _descriptorDict[commandKeys[i]]['type'] == 'num':
                            num_to_create.append(DRP.models.NumMolDescriptorValue(descriptor=descriptorDict[commandKeys[i]], compound=compound, value=resList[i]))
                        elif _descriptorDict[commandKeys[i]]['type'] == 'ord':
                            ord_to_create.append(DRP.models.OrdMolDescriptorValue(descriptor=descriptorDict[commandKeys[i]], compound=compound, value=resList[i]))
                        elif _descriptorDict[commandKeys[i]]['type'] == 'bool':
                            bool_to_create.append(DRP.models.BoolMolDescriptorValue(descriptor=descriptorDict[commandKeys[i]], compound=compound, value=resList[i]))
                        # NOTE: No categorical descriptors are included yet, and since they are more complicated to code I've left it for the moment.
                        # NOTE: Calculation failure values are not included in the documentation, so I've assumed that it doesn't happen, since we have no way of identifying
                        # for it other than for the database to push it out as a part of validation procedures.
                    DRP.models.NumMolDescriptorValue.objects.bulk_create(num_to_create)
                    DRP.models.OrdMolDescriptorValue.objects.bulk_create(ord_to_create)
                    DRP.models.BoolMolDescriptorValue.objects.bulk_create(bool_to_create)
                    if verbose:
                        print "Done with {} ({}/{})".format(compound, i+1, len(compound_dict))
        else:
            raise RuntimeError("cxcalc returned error: {}".fomrat(resErr))
    else:
        raise RuntimeError("cxcalc returned nonzero return code {}".format(calcProc.returncode))

def _calculate(compound, verbose=False):
    notFound = True
    if notFound and (compound.smiles is not None and compound.smiles != ''):
        lecProc = Popen([settings.CHEMAXON_DIR['15.6'] + 'cxcalc', compound.smiles, 'leconformer'], stdout=PIPE, stderr=PIPE, close_fds=True) # lec = lowest energy conformer
        lecProc.wait()
        if lecProc.returncode == 0:
            lec, lecErr = lecProc.communicate()
            notFound = False  
    if notFound and (compound.INCHI is not None and compound.INCHI != ''):
        lecProc = Popen([settings.CHEMAXON_DIR['15.6'] + 'cxcalc', compound.INCHI, 'leconformer'], stdout=PIPE, stderr=PIPE, close_fds=True) # lec = lowest energy conformer
        lecProc.wait()
        if lecProc.returncode == 0:
            lec, lecErr = lecProc.communicate()
            notFound = False
    if not notFound:
        # -N ih means leave off the header row and id column
        calcProc = Popen([settings.CHEMAXON_DIR['15.6'] + 'cxcalc', '-N', 'ih', lec] + [x for x in chain(*(command.split(' ') for command in cxcalcCommands.values()))], stdout=PIPE, stderr=PIPE, close_fds=True) 
        calcProc.wait()
        if calcProc.returncode == 0:
            res, resErr = calcProc.communicate()
            if not resErr:
                resLines = res.split('\n')
                if len(resLines) == 2: #last line is blank
                    resList = resLines[0].split('\t')
                    commandKeys = cxcalcCommands.keys()

                    num_to_create = []
                    ord_to_create = []
                    bool_to_create = []
                    if len(resList) == len(commandKeys):
                        for i in range(len(resList)):
                            if _descriptorDict[commandKeys[i]]['type'] == 'num':
                                num_to_create.append(DRP.models.NumMolDescriptorValue(descriptor=descriptorDict[commandKeys[i]], compound=compound, value=resList[i]))
                            elif _descriptorDict[commandKeys[i]]['type'] == 'ord':
                                ord_to_create.append(DRP.models.OrdMolDescriptorValue(descriptor=descriptorDict[commandKeys[i]], compound=compound, value=resList[i]))
                            elif _descriptorDict[commandKeys[i]]['type'] == 'bool':
                                bool_to_create.append(DRP.models.BoolMolDescriptorValue(descriptor=descriptorDict[commandKeys[i]], compound=compound, value=resList[i]))
                            # NOTE: No categorical descriptors are included yet, and since they are more complicated to code I've left it for the moment.
                            # NOTE: Calculation failure values are not included in the documentation, so I've assumed that it doesn't happen, since we have no way of identifying
                            # for it other than for the database to push it out as a part of validation procedures.
                        
                        if verbose:
                            print "Creating {} numerical, {} ordinal, {} boolean values".format(len(num_to_create), len(ord_to_create), len(bool_to_create))
                        DRP.models.NumMolDescriptorValue.objects.bulk_create(num_to_create)
                        DRP.models.OrdMolDescriptorValue.objects.bulk_create(ord_to_create)
                        DRP.models.BoolMolDescriptorValue.objects.bulk_create(bool_to_create)
                    else:
                        raise RuntimeError("Number of cxcalc commands ({}) does not match number of results ({})".format(len(commandKeys), len(resList)))
            else:
                warnings.warn("cxcalc returned error: {}".format(resErr))
        else:
            warnings.warn("cxcalc exited with nonzero return code {}".format(calcProc.returncode))
    else:
        warnings.warn("Compound not found")
