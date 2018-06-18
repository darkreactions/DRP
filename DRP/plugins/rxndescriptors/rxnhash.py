"""Module for calculating reaction hash descriptor."""
from .utils import setup
import logging

import xxhash

import DRP
from DRP.chemical_data import elements

logger = logging.getLogger(__name__)

calculatorSoftware = 'DRP_xxhash'

_descriptorDict = {
    'rxnSpaceHash1':
        {
            'type': 'cat',
            'name': 'Hash of reaction reactants to partition reaction space',
            'calculatorSoftware': calculatorSoftware,
            'calculatorSoftwareVersion': '0_02_{}'.format(xxhash.VERSION),
            'permittedValues': []
        }
}

descriptorDict = setup(_descriptorDict)


def calculate_many(reaction_set, verbose=False, whitelist=None):
    """Calculate descriptors for this plugin for an entire set of reactions."""
    if verbose:
        logger.info("Creating descriptor dictionary")
    # We're about to use it and leaving it lazy obscures where time is being
    # spent
    descriptorDict.initialise(descriptorDict.descDict)

    for i, reaction in enumerate(reaction_set):
        if verbose:
            logger.info("Calculating {} ({}/{})".format(reaction,
                                                        i + 1, len(reaction_set)))
        _calculate(reaction, descriptorDict,
                   verbose=verbose, whitelist=whitelist)
    
    return descriptorDict

def calculate(reaction, verbose=False, whitelist=None):
    """Calculate the descriptors for this plugin."""
    if verbose:
        logger.info("Creating descriptor dictionary")
    # We're about to use it and leaving it lazy obscures where time is being
    # spent
    descriptorDict.initialise(descriptorDict.descDict)
    _calculate(reaction, descriptorDict, verbose=verbose, whitelist=whitelist)


def _calculate(reaction, descriptorDict, verbose=False, whitelist=None):
    """Calculate descriptors for this plugin with descriptorDict already created."""
    # descriptor Value classes
    cat = DRP.models.CatRxnDescriptorValue
    perm = DRP.models.CategoricalDescriptorPermittedValue

    # reaction space descriptor
    heading = 'rxnSpaceHash1'
    if whitelist is None or heading in whitelist:
        h = xxhash.xxh64()  # generates a hash
        for reactant in reaction.compounds.order_by('pk'):
            h.update('{0:20d}'.format(reactant.pk))
        p = perm.objects.get_or_create(descriptor=descriptorDict[
                                       heading], value=h.hexdigest())[0]
        cat.objects.update_or_create(defaults={
                                     'value': p}, reaction=reaction, descriptor=descriptorDict['rxnSpaceHash1'])[0]
