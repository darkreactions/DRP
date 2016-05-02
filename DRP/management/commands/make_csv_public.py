from DRP.models import PerformedReaction
from django.core.management.base import BaseCommand
import csv
import reimport_reactions
import warning

class Command(BaseCommand):
    def add_arguments(self, parser):
        parser.add_argument('filename', help='Path to csv file.')

    def handle(self, *args, **kwargs):
        fn = kwargs['filename']
        with open(fn) as f:
            reader = csv.reader(f)
            reader.next() # skip headers
        
            for i, row in enumerate(reader):
                ref = reimport_reactions.convert_legacy_reference(row[0])
                self.stdout.write('{}: Making reference {} public'.format(i, ref))
                ps = PerformedReaction.objects.filter(reference=ref)
                if not ps.count() != 1:
                    raise RuntimeError('Found {} reactions with reference {}'.format(ps.count(), ref))
                else:
                    ps.update(public=True)
        
