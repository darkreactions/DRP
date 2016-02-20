
from django.core.management.base import BaseCommand
from DRP.Email import EmailToAdmins
from DRP.models import Reaction
import xxhash


class Command(BaseCommand):

    """Runs a check on the database to make sure reaction hashes don't have collisions.

    The current hash is the same as the one used in rxnSpaceHash1 descriptor.
    Exits nonzero if a collision is found.
    """

    help = 'Checks Reaction Chemical Space Hashes for collisions, which could cause problematic behaviour in model building'

    def handle(self, *args, **kwargs):
        hashDictionary = {}
        collisionCount = 0
        for reaction in Reaction.objects.all():
            reactantString = ''
            h = xxhash.xxh64()
            for reactant in reaction.compounds:
                h.update(reactant.abbrev)
                reactantString += reactant.abbrev
            if h in hashDictionary:
                if hashDictionary[h].hexdigest() != reactantString:
                    collisionCount += 1
            else:
                hashDictionary[h] = reactantString
        if collisionCount > 0:
            e = EmailToAdmins('Dark Reactions Project: Hash Collision Failure',
                              'A collision between reaction space hashes has occured. Please contact the DRP development team and file a bug report.')
            e.send()
            exit(1)
