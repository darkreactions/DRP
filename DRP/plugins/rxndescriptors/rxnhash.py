"""Module for calculating reaction hash descriptor."""
import DRP
from utils import setup
import xxhash

elements = DRP.chemical_data.elements

calculatorSoftware = 'DRP_xxhash'

_descriptorDict = {
    'rxnSpaceHash1':
        {
            'type': 'cat',
            'name': 'Hash of reaction reactants to partition reaction space',
            'calculatorSoftware': calculatorSoftware,
            'calculatorSoftwareVersion': '0.02_{}'.format(xxhash.VERSION),
            'permittedValues': []
        }
}

descriptorDict = setup(_descriptorDict)


def calculate_many(reaction_set, verbose=False):
    """Calculate descriptors for this plugin for an entire set of reactions."""
    if verbose:
        print "Creating descriptor dictionary"
    descriptorDict.initialise(descriptorDict.descDict)  # We're about to use it and leaving it lazy obscures where time is being spent

    for i, reaction in enumerate(reaction_set):
        if verbose:
            print "Calculating {} ({}/{})".format(reaction, i + 1, len(reaction_set))
        _calculate(reaction, descriptorDict, verbose=verbose)


def calculate(reaction, verbose=False):
    """Calculate the descriptors for this plugin."""
    if verbose:
        print "Creating descriptor dictionary"
    descriptorDict.initialise(descriptorDict.descDict)  # We're about to use it and leaving it lazy obscures where time is being spent
    _calculate(reaction, descriptorDict, verbose=verbose)


def _calculate(reaction, descriptorDict, verbose=False):
    """Calculate descriptors for this plugin with descriptorDict already created."""
    # descriptor Value classes
    cat = DRP.models.CatRxnDescriptorValue
    perm = DRP.models.CategoricalDescriptorPermittedValue

    # reaction space descriptor
    h = xxhash.xxh64()  # generates a hash
    for reactant in reaction.compounds.all():
        h.update(reactant.abbrev)
    p = perm.objects.get_or_create(descriptor=descriptorDict['rxnSpaceHash1'], value=h.hexdigest())[0]
    c = cat.objects.get_or_create(reaction=reaction, descriptor=descriptorDict['rxnSpaceHash1'])[0]
    c.value = p
    c.save()
