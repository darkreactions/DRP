from DRP.models import PerformedReaction
from django.core.management.base import BaseCommand
import csv
import reimport_reactions

class Command(BaseCommand):
    help='Make reactions from a csv public. References should be in the first column'
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
                ps = PerformedReaction.objects.filter(convertedLegacyRef=ref)
                if ps.count() == 0:
                    if ref.startswith('xxx'):
                        unmunged_ref = ref + '0'
                        self.stdout.write('{}: Making UNMUNGED reference {} public'.format(i, unmunged_ref))
                        ps = PerformedReaction.objects.filter(convertedLegacyRef=unmunged_ref)
                        if ps.count() != 1:
                            ps.update(public=True)
                    else:
                        raise RuntimeError('Found {} reactions with reference {}'.format(ps.count(), ref))
                else:
                    ps.update(public=True)
        
