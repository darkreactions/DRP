"""Used for deleting reactions specified in a tsv."""
from django.core.management.base import BaseCommand
from django.db import transaction
from os import path
import csv
from DRP.models import PerformedReaction
import reimport_reactions


class Command(BaseCommand):
    
    """Used for deleting reactions specified in a tsv."""
    
    help = 'Deletes reactions from a tsv file'

    def add_arguments(self, parser):
        """Add arguments ot the parser."""
        parser.add_argument('directory', help='The directory where the tsv files are')

    def handle(self, *args, **kwargs):
        """Handle the command call."""
        folder = kwargs['directory']
        self.stdout.write('Deleting reactions')

        with transaction.atomic():
            with open(path.join(folder, 'performedReactions.tsv')) as reactions:
                reader = csv.DictReader(reactions, delimiter='\t')
                for i, r in enumerate(reader):
                    ref = r['reference']
                    if ref.lower() != reimport_reactions.convert_legacy_reference(ref):
                        ps = PerformedReaction.objects.filter(reference=ref.lower())
                        if ps:
                            self.stdout.write('{}: Deleting reaction with converted legacy reference {}'.format(i, ref))
                            ps.delete()
