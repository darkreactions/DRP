"""Basic reaction descriptors calculation module."""
import logging
from itertools import chain

import xxhash
from numpy import mean, average as wmean
from scipy.stats import gmean
from django.db.models import Sum

import DRP
from DRP.plugins.moldescriptors.chemaxon import _pHDependentDescriptors
from .utils import setup
from DRP.chemical_data import elements

logger = logging.getLogger(__name__)

calculatorSoftware = 'DRP'
# number of values to create at a time. Should probably be <= 5000
create_threshold = 5000

_descriptorDict = {}

# The following adds descriptors to the dictionary in an automated way to
# save on voluminous code

for element in elements:
    _descriptorDict[element + '_mols'] = {
        'type': 'num',
        'name': 'Mols of {} in a reaction.'.format(element),
        'calculatorSoftware': calculatorSoftware,
        'calculatorSoftwareVersion': '0_02',
        'maximum': None,
        'minimum': 0,
    }


# descriptors for generalised aggregation across compound roles
def make_dict():
    """Make a dictionary of descriptors based on known patterns."""
    _reaction_pH_Descriptors = {}
    weightings = ('molarity', 'count')
    for compoundRole in DRP.models.CompoundRole.objects.all():
        for w in weightings:
            _descriptorDict['{}_amount_{}'.format(compoundRole.label, w)] = {
                'type': 'num',
                'name': 'Amount of compounds in the reaction belonging in the role "{}" weighted by {}'.format(compoundRole.label, w),
                'calculatorSoftware': calculatorSoftware,
                'calculatorSoftwareVersion': '0_02',
                'maximum': None,
                'minimum': 0,
            }
        for descriptor in DRP.models.CatMolDescriptor.objects.all():
            for w in weightings:
                for permValue in descriptor.permittedValues.all():
                    _descriptorDict['{}_{}_{}_{}'.format(compoundRole.label, descriptor.csvHeader, permValue.value, w)] = {
                        'type': 'num',
                        'name': 'Amount of reactants in category {} for descriptor "{}" weighted by reactant {} in compound role {}.'.format(
                                permValue.value, descriptor.name, w, compoundRole.label),
                        'calculatorSoftware': calculatorSoftware,
                        'calculatorSoftwareVersion': '0_02',
                        'maximum': None,
                        'minimum': 0,
                    }
        for descriptor in DRP.models.OrdMolDescriptor.objects.all():
            for w in weightings:
                # because python...
                for i in range(descriptor.minimum, descriptor.maximum + 1):
                    _descriptorDict['{}_{}_{}_{}'.format(compoundRole.label, descriptor.csvHeader, i, w)] = {
                        'type': 'num',
                        'name': 'Amount of reactants with value {} for descriptor "{}" weighted by reactant {} in compound role {}.'.format(
                            i, descriptor.name, w, compoundRole.label),
                        'calculatorSoftware': calculatorSoftware,
                        'calculatorSoftwareVersion': '0_02',
                        'maximum': None,
                        'minimum': 0,
                    }
        for descriptor in DRP.models.BoolMolDescriptor.objects.all():
            for w in weightings:
                for value in ('True', 'False'):
                    _descriptorDict['{}_{}_{}_{}'.format(compoundRole.label, descriptor.csvHeader, value, w)] = {
                        'type': 'num',
                        'name': 'Amount of reactants with value {} for descriptor "{}" weighted by reactant {} in compound role {}.'.format(
                            value, descriptor.name, w, compoundRole.label),
                        'calculatorSoftware': calculatorSoftware,
                        'calculatorSoftwareVersion': '0_02',
                        'maximum': None,
                        'minimum': 0,
                    }

            _descriptorDict['{}_{}_any'.format(compoundRole.label, descriptor.csvHeader)] = {
                'type': 'bool',
                'name': 'Whether any reactants have value True for descriptor "{}" in compound role {}.'.format(
                        descriptor.name, compoundRole.label),
                'calculatorSoftware': calculatorSoftware,
                'calculatorSoftwareVersion': '1_5',
            }

        for descriptor in DRP.models.NumMolDescriptor.objects.all():
            _descriptorDict['{}_{}_{}'.format(compoundRole.label, descriptor.csvHeader, 'Max')] = {
                'type': 'num',
                'name': 'Maximum value for {} aggregated across compounds in role "{}"'.format(descriptor.name, compoundRole.label),
                'calculatorSoftware': calculatorSoftware,
                'calculatorSoftwareVersion': '0_02',
                'maximum': descriptor.maximum,
                'minimum': descriptor.minimum,
            }
            _descriptorDict['{}_{}_{}'.format(compoundRole.label, descriptor.csvHeader, 'Range')] = {
                'type': 'num',
                'name': 'Range for {} aggregated across compounds in role "{}"'.format(descriptor.name, compoundRole.label),
                'calculatorSoftware': calculatorSoftware,
                'calculatorSoftwareVersion': '0_02',
                'maximum': (descriptor.maximum - descriptor.minimum) if (descriptor.maximum is not None and descriptor.minimum is not None) else None,
                'minimum': 0,
            }
            for w in weightings:
                _descriptorDict['{}_{}_{}_{}'.format(compoundRole.label, descriptor.csvHeader, 'gmean', w)] = {
                    'type': 'num',
                    'name': 'Geometric Mean for {} aggregated across compounds in role "{}" normalised by {}'.format(descriptor.name, compoundRole.label, w),
                    'calculatorSoftware': calculatorSoftware,
                    'calculatorSoftwareVersion': '0_02',
                    'maximum': None,
                    'minimum': 0,
                }

        for heading, d in _pHDependentDescriptors.items():
            _reaction_pH_Descriptors['{}_{}_pHreaction_{}_{}_{}'.format(compoundRole.label, heading, d['calculatorSoftware'], d['calculatorSoftwareVersion'], 'Max')] = {
                'type': 'num',
                'name': 'Maximum value for {} aggregated across compounds in role "{}"'.format(d['name'] + ' at reaction pH', compoundRole.label),
                'calculatorSoftware': calculatorSoftware,
                'calculatorSoftwareVersion': '1_5',
                'maximum': d['maximum'],
                'minimum': d['minimum'],
            }
            _reaction_pH_Descriptors['{}_{}_pHreaction_{}_{}_{}'.format(compoundRole.label, heading, d['calculatorSoftware'], d['calculatorSoftwareVersion'], 'Range')] = {
                'type': 'num',
                'name': 'Range for {} aggregated across compounds in role "{}"'.format(d['name'] + ' at reaction pH', compoundRole.label),
                'calculatorSoftware': calculatorSoftware,
                'calculatorSoftwareVersion': '1_5',
                'maximum': d['maximum'] - d['minimum'] if (d['maximum'] is not None and d['minimum'] is not None) else None,
                'minimum': 0,
            }
            for w in weightings:
                _reaction_pH_Descriptors['{}_{}_pHreaction_{}_{}_{}_{}'.format(compoundRole.label, heading, d['calculatorSoftware'], d['calculatorSoftwareVersion'], 'gmean', w)] = {
                    'type': 'num',
                    'name': 'Geometric Mean for {} aggregated across compounds in role "{}" normalised by {}'.format(d['name'] + ' at reaction pH', compoundRole.label, w),
                    'calculatorSoftware': calculatorSoftware,
                    'calculatorSoftwareVersion': '1_5',
                    'maximum': None,
                    'minimum': 0,
                }

    _descriptorDict.update(_reaction_pH_Descriptors)

    descriptorDict = setup(_descriptorDict)
    return descriptorDict, _reaction_pH_Descriptors

# There's a lot of DRY violation here because I was playing with a few different methods.
# We should decide which method we want for deletion and work on
# variations of that.


def delete_descriptors_many(reaction_set, descriptorDict):
    """Bulk deletion of descriptors."""
    # This could be the same as delete descriptors if we're ok with deleting the descriptor even if not
    # descriptorValues.count() == roleQuantities.count() and not any(descriptorValue.value is None for descriptorValue in descriptorValues)
    # I think yes, but wasn't completely sure, so made the separate function

    NumRxnDescriptor = DRP.models.NumRxnDescriptor
    BoolRxnDescriptor = DRP.models.BoolRxnDescriptor
    OrdRxnDescriptor = DRP.models.OrdRxnDescriptor
    CatRxnDescriptor = DRP.models.CatRxnDescriptor

    descriptors_to_delete = []
    for element in elements:
        descriptors_to_delete.append(descriptorDict[element + '_mols'])

    allCompoundQuantities = DRP.models.CompoundQuantity.objects.filter(
        reaction__in=reaction_set)
    for compoundRole in DRP.models.CompoundRole.objects.all():
        roleQuantities = allCompoundQuantities.filter(role=compoundRole)
        descriptors_to_delete.append(
            descriptorDict['{}_amount_count'.format(compoundRole.label)])
        #  number of species in reaction with this role
        descriptors_to_delete.append(
            descriptorDict['{}_amount_molarity'.format(compoundRole.label)])

        if roleQuantities.exists():
            descriptorValues = DRP.models.NumMolDescriptorValue.objects.filter(
                compound__in=[quantity.compound for quantity in roleQuantities])
            for descriptor in DRP.models.NumMolDescriptor.objects.all():
                descriptors_to_delete.append(descriptorDict['{}_{}_{}'.format(
                    compoundRole.label, descriptor.csvHeader, 'Max')])
                descriptors_to_delete.append(descriptorDict['{}_{}_{}'.format(
                    compoundRole.label, descriptor.csvHeader, 'Range')])
                descriptors_to_delete.append(descriptorDict['{}_{}_{}_{}'.format(
                    compoundRole.label, descriptor.csvHeader, 'gmean', 'molarity')])
                descriptors_to_delete.append(descriptorDict['{}_{}_{}_{}'.format(
                    compoundRole.label, descriptor.csvHeader, 'gmean', 'count')])
            for descriptor in DRP.models.OrdMolDescriptor.objects.all():
                # because Still python...
                for i in range(descriptor.minimum, descriptor.maximum + 1):
                    descriptors_to_delete.append(descriptorDict['{}_{}_{}_count'.format(
                        compoundRole.label, descriptor.csvHeader, i)])
                    descriptors_to_delete.append(descriptorDict['{}_{}_{}_molarity'.format(
                        compoundRole.label, descriptor.csvHeader, i)])
            for descriptor in DRP.models.BoolMolDescriptor.objects.all():
                for i in (True, False):  # because Still python...
                    descriptors_to_delete.append(descriptorDict['{}_{}_{}_count'.format(
                        compoundRole.label, descriptor.csvHeader, i)])
                    descriptors_to_delete.append(descriptorDict['{}_{}_{}_molarity'.format(
                        compoundRole.label, descriptor.csvHeader, i)])
                descriptors_to_delete.append(
                    descriptorDict['{}_{}_any'.format(compoundRole.label, descriptor.csvHeader)])
            for descriptor in DRP.models.CatMolDescriptor.objects.all():
                for permValue in descriptor.permittedValues.all():
                    descriptors_to_delete.append(descriptorDict['{}_{}_{}_count'.format(
                        compoundRole.label, descriptor.csvHeader, permValue.value)])
                    descriptors_to_delete.append(descriptorDict['{}_{}_{}_molarity'.format(
                        compoundRole.label, descriptor.csvHeader, permValue.value)])

    _delete_values(reaction_set, descriptors_to_delete)


def delete_descriptors(reaction, descriptorDict, whitelist=None):
    """Deletion of descriptors with some restrictions."""
    allCompoundQuantities = DRP.models.CompoundQuantity.objects.filter(
        reaction=reaction)
    NumRxnDescriptor = DRP.models.NumRxnDescriptor
    BoolRxnDescriptor = DRP.models.BoolRxnDescriptor
    OrdRxnDescriptor = DRP.models.OrdRxnDescriptor
    CatRxnDescriptor = DRP.models.CatRxnDescriptor

    descriptors_to_delete = []
    for element in elements:
        descriptors_to_delete.append(descriptorDict[element + '_mols'])

    for compoundRole in DRP.models.CompoundRole.objects.all():
        roleQuantities = allCompoundQuantities.filter(role=compoundRole)
        descriptors_to_delete.append(
            descriptorDict['{}_amount_count'.format(compoundRole.label)])
        #  number of species in reaction with this role
        descriptors_to_delete.append(
            descriptorDict['{}_amount_molarity'.format(compoundRole.label)])

        if roleQuantities.exists():
            for descriptor in DRP.models.NumMolDescriptor.objects.all():
                descriptorValues = DRP.models.NumMolDescriptorValue.objects.filter(
                    compound__in=[quantity.compound for quantity in roleQuantities], descriptor=descriptor)
                # Only do the calculation if the right number of descriptor
                # values are present and all of them are not NULL
                if descriptorValues.count() == roleQuantities.count() and not any(descriptorValue.value is None for descriptorValue in descriptorValues):
                    descriptors_to_delete.append(descriptorDict['{}_{}_{}'.format(
                        compoundRole.label, descriptor.csvHeader, 'Max')])
                    descriptors_to_delete.append(descriptorDict['{}_{}_{}'.format(
                        compoundRole.label, descriptor.csvHeader, 'Range')])
                    descriptors_to_delete.append(descriptorDict['{}_{}_{}_{}'.format(
                        compoundRole.label, descriptor.csvHeader, 'gmean', 'molarity')])
                    descriptors_to_delete.append(descriptorDict['{}_{}_{}_{}'.format(
                        compoundRole.label, descriptor.csvHeader, 'gmean', 'count')])
            for descriptor in DRP.models.OrdMolDescriptor.objects.all():
                descriptorValues = DRP.models.OrdMolDescriptorValue.objects.filter(
                    compound__in=[quantity.compound for quantity in roleQuantities], descriptor=descriptor)
                # Only do the calculation if the right number of descriptor
                # values are present and all of them are not NULL
                if descriptorValues.count() == roleQuantities.count() and not any(descriptorValue.value is None for descriptorValue in descriptorValues):
                    # because Still python...
                    for i in range(descriptor.minimum, descriptor.maximum + 1):
                        descriptors_to_delete.append(descriptorDict['{}_{}_{}_count'.format(
                            compoundRole.label, descriptor.csvHeader, i)])
                        descriptors_to_delete.append(descriptorDict['{}_{}_{}_molarity'.format(
                            compoundRole.label, descriptor.csvHeader, i)])
            for descriptor in DRP.models.BoolMolDescriptor.objects.all():
                descriptorValues = DRP.models.BoolMolDescriptorValue.objects.filter(
                    compound__in=[quantity.compound for quantity in roleQuantities], descriptor=descriptor)
                if descriptorValues.count() == roleQuantities.count() and not any(descriptorValue.value is None for descriptorValue in descriptorValues):
                    for i in (True, False):  # because Still python...
                        descriptors_to_delete.append(descriptorDict['{}_{}_{}_count'.format(
                            compoundRole.label, descriptor.csvHeader, i)])
                        descriptors_to_delete.append(descriptorDict['{}_{}_{}_molarity'.format(
                            compoundRole.label, descriptor.csvHeader, i)])
                    descriptors_to_delete.append(
                        descriptorDict['{}_{}_any'.format(compoundRole.label, descriptor.csvHeader)])
            for descriptor in DRP.models.CatMolDescriptor.objects.all():
                descriptorValues = DRP.models.CatMolDescriptorValue.objects.filter(
                    compound__in=[quantity.compound for quantity in roleQuantities], descriptor=descriptor)
                if descriptorValues.count() == roleQuantities.count() and not any(descriptorValue.value is None for descriptorValue in descriptorValues):
                    for permValue in descriptor.permittedValues.all():
                        descriptors_to_delete.append(descriptorDict['{}_{}_{}_count'.format(
                            compoundRole.label, descriptor.csvHeader, permValue.value)])
                        descriptors_to_delete.append(descriptorDict['{}_{}_{}_molarity'.format(
                            compoundRole.label, descriptor.csvHeader, permValue.value)])

    if whitelist is not None:
        descriptors_to_delete = [
            desc for desc in descriptors_to_delete if desc in whitelist]
    _delete_values([reaction], descriptors_to_delete)


def _delete_values(reaction_set, descriptors_to_delete):
    NumRxnDescriptor = DRP.models.NumRxnDescriptor
    BoolRxnDescriptor = DRP.models.BoolRxnDescriptor
    OrdRxnDescriptor = DRP.models.OrdRxnDescriptor
    CatRxnDescriptor = DRP.models.CatRxnDescriptor

    if len(reaction_set) == 1:
        rxn = reaction_set[0]
        DRP.models.NumRxnDescriptorValue.objects.filter(reaction=rxn, descriptor__in=[
                                                        desc for desc in descriptors_to_delete if isinstance(desc, NumRxnDescriptor)]).delete()
        DRP.models.BoolRxnDescriptorValue.objects.filter(reaction=rxn, descriptor__in=[
                                                         desc for desc in descriptors_to_delete if isinstance(desc, BoolRxnDescriptor)]).delete()
        DRP.models.OrdRxnDescriptorValue.objects.filter(reaction=rxn, descriptor__in=[
                                                        desc for desc in descriptors_to_delete if isinstance(desc, OrdRxnDescriptor)]).delete()
        DRP.models.CatRxnDescriptorValue.objects.filter(reaction=rxn, descriptor__in=[
                                                        desc for desc in descriptors_to_delete if isinstance(desc, CatRxnDescriptor)]).delete()
    else:
        DRP.models.NumRxnDescriptorValue.objects.filter(reaction__in=reaction_set, descriptor__in=[
                                                        desc for desc in descriptors_to_delete if isinstance(desc, NumRxnDescriptor)]).delete()
        DRP.models.BoolRxnDescriptorValue.objects.filter(reaction__in=reaction_set, descriptor__in=[
                                                         desc for desc in descriptors_to_delete if isinstance(desc, BoolRxnDescriptor)]).delete()
        DRP.models.OrdRxnDescriptorValue.objects.filter(reaction__in=reaction_set, descriptor__in=[
                                                        desc for desc in descriptors_to_delete if isinstance(desc, OrdRxnDescriptor)]).delete()
        DRP.models.CatRxnDescriptorValue.objects.filter(reaction__in=reaction_set, descriptor__in=[
                                                        desc for desc in descriptors_to_delete if isinstance(desc, CatRxnDescriptor)]).delete()


def calculate_many(reaction_set, verbose=False, bulk_delete=False, whitelist=None):
    """Calculate descriptors for this plugin for an entire set of reactions."""
    if verbose:
        logger.info("Creating descriptor dictionary")
    descriptorDict, _reaction_pH_Descriptors = make_dict()
    # We're about to use it and leaving it lazy obscures where time is being
    # spent
    descriptorDict.initialise(descriptorDict.descDict)

    if whitelist is None:
        descs_to_delete = descriptorDict.values()
    else:
        descs_to_delete = [descriptorDict[k]
                           for k in descriptorDict.keys() if k in whitelist]

    if bulk_delete:
        if verbose:
            logger.info("Deleting all old descriptor values")
        _delete_values(reaction_set, descs_to_delete)

    num_vals_to_create = []
    bool_vals_to_create = []
    for i, reaction in enumerate(reaction_set):
        if verbose:
            logger.info("{} ({}/{})".format(reaction, i + 1, len(reaction_set)))
        if not bulk_delete:
            if verbose:
                logger.info("Deleting old descriptor values")
            _delete_values([reaction], descs_to_delete)

        if verbose:
            logger.info("Calculating new values.".format(reaction, i + 1, len(reaction_set)))
        num_vals_to_create, bool_vals_to_create = _calculate(
            reaction, descriptorDict, verbose=verbose, whitelist=whitelist, num_vals_to_create=num_vals_to_create, bool_vals_to_create=bool_vals_to_create)

        if len(num_vals_to_create) > create_threshold:
            if verbose:
                logger.info("Creating {} Numeric values".format(len(num_vals_to_create)))
            DRP.models.NumRxnDescriptorValue.objects.bulk_create(
                num_vals_to_create)
            num_vals_to_create = []
        if len(bool_vals_to_create) > create_threshold:
            if verbose:
                logger.info("Creating {} Boolean values".format(len(bool_vals_to_create)))
            DRP.models.BoolRxnDescriptorValue.objects.bulk_create(
                bool_vals_to_create)
            bool_vals_to_create = []

    if verbose:
        logger.info("Creating {} Numeric values".format(len(num_vals_to_create)))
    DRP.models.NumRxnDescriptorValue.objects.bulk_create(num_vals_to_create)
    if verbose:
        logger.info("Creating {} Boolean values".format(len(bool_vals_to_create)))
    DRP.models.BoolRxnDescriptorValue.objects.bulk_create(bool_vals_to_create)

    if verbose:
        logger.info("Creating reaction pH values")

    num_vals_to_create = []
    for i, reaction in enumerate(reaction_set):
        if verbose:
            logger.info("{} ({}/{})".format(reaction, i + 1, len(reaction_set)))

        if verbose:
            logger.info("Calculating new values.".format(reaction, i + 1, len(reaction_set)))
        num_vals_to_create = _calculateRxnpH(reaction, descriptorDict, _reaction_pH_Descriptors,
                                             verbose=verbose, whitelist=whitelist, vals_to_create=num_vals_to_create)

        if len(num_vals_to_create) > create_threshold:
            if verbose:
                logger.info("Creating {} Numeric values".format(len(num_vals_to_create)))
            DRP.models.NumRxnDescriptorValue.objects.bulk_create(
                num_vals_to_create)
            num_vals_to_create = []

    if verbose:
        logger.info("Creating {} Numeric values".format(len(num_vals_to_create)))
    DRP.models.NumRxnDescriptorValue.objects.bulk_create(num_vals_to_create)


def calculate(reaction, verbose=False, whitelist=None):
    """Calculate the descriptors for this plugin."""
    # Set up the actual descriptor dictionary.
    if verbose:
        logger.info("Creating descriptor dictionary")
    descriptorDict = make_dict()[0]
    if whitelist is None:
        descs_to_delete = descriptorDict.values()
    else:
        descs_to_delete = [k[v]
                           for k in descriptorDict.keys() if k in whitelist]
    if verbose:
        logger.info("Deleting old values")
    _delete_values([reaction], descs_to_delete)
    if verbose:
        logger.info("Calculating new values")
    num_vals_to_create, bool_vals_to_create = _calculate(
        reaction, descriptorDict, verbose=verbose, whitelist=whitelist)

    if verbose:
        logger.info("Creating {} Numeric values".format(len(num_vals_to_create)))
    DRP.models.NumRxnDescriptorValue.objects.bulk_create(num_vals_to_create)
    if verbose:
        logger.info("Creating {} Boolean values".format(len(bool_vals_to_create)))
    DRP.models.BoolRxnDescriptorValue.objects.bulk_create(bool_vals_to_create)

    if verbose:
        logger.info("Calculating reaction pH values")
    num_vals_to_create = _calculateRxnpH(reaction, descriptorDict, _reaction_pH_descriptors,
                                         verbose=verbose, whitelist=whitelist, vals_to_create=num_vals_to_create)

    if verbose:
        logger.info("Creating {} Numeric values".format(len(num_vals_to_create)))
    DRP.models.NumRxnDescriptorValue.objects.bulk_create(num_vals_to_create)


def _calculate(reaction, descriptorDict, verbose=False, whitelist=None, num_vals_to_create=[], bool_vals_to_create=[]):
    """Calculate with the descriptorDict already created and previous descriptor values deleted."""
    # descriptor Value classes
    CompoundQuantity = DRP.models.CompoundQuantity
    num = DRP.models.NumRxnDescriptorValue
    cat = DRP.models.CatRxnDescriptorValue
    perm = DRP.models.CategoricalDescriptorPermittedValue

    # Calculate the elemental molarities
    allCompoundQuantities = CompoundQuantity.objects.filter(reaction=reaction)

    for element in elements:
        heading = element + '_mols'
        if whitelist is None or heading in whitelist:
            n = num(
                reaction=reaction,
                descriptor=descriptorDict[heading],
            )
            if any(quantity.amount is None for quantity in allCompoundQuantities):
                n.value = None
            else:
                n.value = float(sum((float(quantity.compound.elements[element]['stoichiometry']) * float(
                    quantity.amount) if element in quantity.compound.elements else 0) for quantity in allCompoundQuantities))

            num_vals_to_create.append(n)

    for compoundRole in DRP.models.CompoundRole.objects.all():
        roleQuantities = allCompoundQuantities.filter(role=compoundRole)
        heading = '{}_amount_count'.format(compoundRole.label)
        if whitelist is None or heading in whitelist:
            #  number of species in reaction with this role
            n = num(
                reaction=reaction,
                descriptor=descriptorDict[heading],
            )
            n.value = roleQuantities.count()
            num_vals_to_create.append(n)

        #  moles of reactant filling this role in this reaction

        if any(quantity.amount is None for quantity in roleQuantities):
            roleMoles = None
        else:
            roleMoles = sum(quantity.amount for quantity in roleQuantities)
        heading = '{}_amount_molarity'.format(compoundRole.label)
        if whitelist is None or heading in whitelist:
            n = num(
                reaction=reaction,
                descriptor=descriptorDict[heading],
            )
            n.value = roleMoles
            num_vals_to_create.append(n)

        if roleQuantities.exists():
            for descriptor in DRP.models.NumMolDescriptor.objects.all():
                descriptorValues = DRP.models.NumMolDescriptorValue.objects.filter(
                    compound__in=[quantity.compound for quantity in roleQuantities], descriptor=descriptor)
                # Only do the calculation if the right number of descriptor values are present and all of them are not NULL
                # this count now takes longer than actually creating the value. Can we remove it? (Same below)
                # Ah looks like this is silently skipping inorganic properties for organics. That seems like a bad way to do it...
                # Means I have to silence the warnings below (which previously lead me to a whole bunch of uncalculated stuff.
                # I really don't like failing silently and it seems to have
                # gotten us into some trouble

                if descriptorValues.count() == roleQuantities.count() and not any(descriptorValue.value is None for descriptorValue in descriptorValues):
                    heading = '{}_{}_{}'.format(
                        compoundRole.label, descriptor.csvHeader, 'Max')
                    if whitelist is None or heading in whitelist:
                        n = num(
                            reaction=reaction,
                            descriptor=descriptorDict[heading],
                        )
                        if any(descriptorValue.value is None for descriptorValue in descriptorValues):
                            n.value = None
                        else:
                            n.value = max(
                                descriptorValue.value for descriptorValue in descriptorValues)
                        num_vals_to_create.append(n)

                    heading = '{}_{}_{}'.format(
                        compoundRole.label, descriptor.csvHeader, 'Range')
                    if whitelist is None or heading in whitelist:
                        n = num(
                            reaction=reaction,
                            descriptor=descriptorDict[heading],
                        )
                        if any(descriptorValue.value is None for descriptorValue in descriptorValues):
                            n.value = None
                        else:
                            n.value = max(descriptorValue.value for descriptorValue in descriptorValues) - min(
                                descriptorValue.value for descriptorValue in descriptorValues)
                        num_vals_to_create.append(n)

                    heading = '{}_{}_{}_{}'.format(
                        compoundRole.label, descriptor.csvHeader, 'gmean', 'molarity')
                    if whitelist is None or heading in whitelist:
                        n = num(
                            reaction=reaction,
                            descriptor=descriptorDict[heading],
                        )
                        if any(descriptorValues.get(compound=quantity.compound).value is None for quantity in roleQuantities) or roleMoles == 0 or roleMoles is None:
                            n.value = None
                        elif any(descriptorValues.get(compound=quantity.compound).value == 0 for quantity in roleQuantities):
                            n.value = 0
                        elif any(descriptorValues.get(compound=quantity.compound).value < 0 for quantity in roleQuantities):
                            raise ValueError(
                                'Cannot take geometric mean of negative values. This descriptor ({}) should not use a geometric mean.'.format(descriptor))
                        else:
                            n.value = gmean(list(descriptorValues.get(compound=quantity.compound).value * float(
                                quantity.amount / roleMoles) for quantity in roleQuantities))
                        num_vals_to_create.append(n)

                    heading = '{}_{}_{}_{}'.format(
                        compoundRole.label, descriptor.csvHeader, 'gmean', 'count')
                    if whitelist is None or heading in whitelist:
                        n = num(
                            reaction=reaction,
                            descriptor=descriptorDict[heading],
                        )
                        if any(descriptorValues.get(compound=quantity.compound).value is None for quantity in roleQuantities) or roleMoles == 0 or roleMoles is None:
                            n.value = None
                        elif any(descriptorValues.get(compound=quantity.compound).value == 0 for quantity in roleQuantities):
                            n.value = 0
                        elif any(descriptorValues.get(compound=quantity.compound).value < 0 for quantity in roleQuantities):
                            raise ValueError(
                                'Cannot take geometric mean of negative values. This descriptor ({}) should not use a geometric mean.'.format(descriptor))
                        else:
                            n.value = gmean(list(descriptorValues.get(
                                compound=quantity.compound).value for quantity in roleQuantities))
                        num_vals_to_create.append(n)
                # elif descriptorValues.count() != roleQuantities.count():
                    # logger.warning("Skipping {} because there are {} descriptorValues and {} roleQuantities".format(descriptor.heading, descriptorValues.count(), roleQuantities.count()))
                # else:
                    # logger.warning("Skipping {} because some descriptorValues are None".format(descriptor.heading))
            for descriptor in DRP.models.OrdMolDescriptor.objects.all():
                #  Only do the calculation if the right number of descriptor values are present and all of them are not NULL
                # this count now takes longer than actually creating
                # the value. Can we remove it? (Same below)
                descriptorValues = DRP.models.OrdMolDescriptorValue.objects.filter(
                    compound__in=[quantity.compound for quantity in roleQuantities], descriptor=descriptor)
                if descriptorValues.count() == roleQuantities.count() and not any(descriptorValue.value is None for descriptorValue in descriptorValues):
                    # because Still python...
                    for i in range(descriptor.minimum, descriptor.maximum + 1):
                        heading = '{}_{}_{}_count'.format(
                            compoundRole.label, descriptor.csvHeader, i)
                        if whitelist is None or heading in whitelist:
                            n = num(
                                reaction=reaction,
                                descriptor=descriptorDict[heading],
                            )
                            n.value = sum(
                                1 for value in descriptorValues if value.value == i)
                            num_vals_to_create.append(n)

                        heading = '{}_{}_{}_molarity'.format(
                            compoundRole.label, descriptor.csvHeader, i)
                        if whitelist is None or heading in whitelist:
                            n = num(
                                reaction=reaction,
                                descriptor=descriptorDict[heading],
                            )
                            quantities = roleQuantities.filter(
                                compound__ordmoldescriptorvalue__value=i, compound__ordmoldescriptorvalue__descriptor__pk=descriptor.pk)
                            if any(quantity.amount is None for quantity in quantities):
                                n.value = None
                            else:
                                n.value = sum(
                                    quantity.amount for quantity in quantities)
                            num_vals_to_create.append(n)
            for descriptor in DRP.models.BoolMolDescriptor.objects.all():
                descriptorValues = DRP.models.BoolMolDescriptorValue.objects.filter(
                    compound__in=[quantity.compound for quantity in roleQuantities], descriptor=descriptor)
                if descriptorValues.count() == roleQuantities.count() and not any(descriptorValue.value is None for descriptorValue in descriptorValues):
                    for i in (True, False):  # because Still python...
                        heading = '{}_{}_{}_count'.format(
                            compoundRole.label, descriptor.csvHeader, i)
                        if whitelist is None or heading in whitelist:
                            n = num(
                                reaction=reaction,
                                descriptor=descriptorDict[heading],
                            )
                            n.value = sum(
                                1 for value in descriptorValues if value.value == i)
                            num_vals_to_create.append(n)

                        heading = '{}_{}_{}_molarity'.format(
                            compoundRole.label, descriptor.csvHeader, i)
                        if whitelist is None or heading in whitelist:
                            n = num(
                                reaction=reaction,
                                descriptor=descriptorDict[heading],
                            )
                            quantities = roleQuantities.filter(
                                compound__boolmoldescriptorvalue__value=i, compound__boolmoldescriptorvalue__descriptor__pk=descriptor.pk)
                            if any(quantity.amount is None for quantity in quantities):
                                n.value = None
                            else:
                                n.value = sum(
                                    quantity.amount for quantity in quantities)
                            num_vals_to_create.append(n)
                    heading = '{}_{}_any'.format(
                        compoundRole.label, descriptor.csvHeader)
                    if whitelist is None or heading in whitelist:
                        b = DRP.models.BoolRxnDescriptorValue(
                            reaction=reaction,
                            descriptor=descriptorDict[heading],
                            value=any(descriptorValues)
                        )
                        bool_vals_to_create.append(b)

            for descriptor in DRP.models.CatMolDescriptor.objects.all():
                descriptorValues = DRP.models.CatMolDescriptorValue.objects.filter(
                    compound__in=[quantity.compound for quantity in roleQuantities], descriptor=descriptor)
                if descriptorValues.count() == roleQuantities.count() and not any(descriptorValue.value is None for descriptorValue in descriptorValues):
                    for permValue in descriptor.permittedValues.all():
                        heading = '{}_{}_{}_count'.format(
                            compoundRole.label, descriptor.csvHeader, permValue.value)
                        if whitelist is None or heading in whitelist:
                            n = num(
                                reaction=reaction,
                                descriptor=descriptorDict[heading],
                            )
                            n.value = descriptorValues.filter(
                                value=permValue).count()
                            num_vals_to_create.append(n)
                        heading = '{}_{}_{}_molarity'.format(
                            compoundRole.label, descriptor.csvHeader, permValue.value)
                        if whitelist is None or heading in whitelist:
                            n = num(
                                reaction=reaction,
                                descriptor=descriptorDict[heading],
                            )
                            quantities = roleQuantities.filter(
                                compound__catmoldescriptorvalue__value=permValue, compound__catmoldescriptorvalue__descriptor__pk=descriptor.pk)
                            if any(quantity.amount is None for quantity in quantities):
                                n.value = None
                            else:
                                n.value = sum(
                                    quantity.amount for quantity in quantities)
                            num_vals_to_create.append(n)
    return num_vals_to_create, bool_vals_to_create


def _calculateRxnpH(reaction, descriptorDict, _reaction_pH_Descriptors, verbose=False, whitelist=None, vals_to_create=[]):
    try:
        reaction_pH = DRP.models.NumRxnDescriptorValue.objects.get(
            reaction=reaction, descriptor__heading='reaction_pH').value
    except DRP.models.NumRxnDescriptorValue.DoesNotExist:
        logger.warning(
            'Reaction {} has no pH value. Cannot create reaction pH descriptors'.format(reaction))
    else:
        reaction_pH_string = str(reaction_pH).replace(
            '.', '_')  # R compatibility
        for heading in _reaction_pH_Descriptors.keys():
            d = descriptorDict[heading]
            if reaction_pH is None:
                pH_descriptor_value = DRP.models.NumRxnDescriptorValue(
                    descriptor=d, reaction=reaction, value=None)
            else:
                reaction_pH_descriptor_heading = d.heading.replace(
                    '_pHreaction_', '_pH{}_'.format(reaction_pH_string))
                try:
                    pH_descriptor_value = DRP.models.NumRxnDescriptorValue.objects.get(
                        descriptor__heading=reaction_pH_descriptor_heading, reaction=reaction).value
                except DRP.models.NumRxnDescriptorValue.DoesNotExist:
                    # this creates a ton of warnings about pH and Ox roles because those mostly don't exist
                    # Not sure that silence wouldn't be better, but I've been burned too many times by silent failure
                    # Came up with this compromise. The warnings are identical
                    # so the default python settings squelch them
                    if heading.startswith('pH_') or heading.startswith('Ox_'):
                        logger.warning(
                            'Could not find descriptor value for a pH or Ox role descriptor.')
                    else:
                        logger.warning('Could not find descriptor value for descriptor {} and reaction {}'.format(
                            reaction_pH_descriptor_heading, reaction))
                else:
                    reaction_pH_dv = DRP.models.NumRxnDescriptorValue(
                        descriptor=d, reaction=reaction, value=pH_descriptor_value)
                    vals_to_create.append(reaction_pH_dv)

    return vals_to_create