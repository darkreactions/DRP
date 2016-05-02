from DRP.models import PerformedReaction
from django.core.management.base import BaseCommand

PerformedReaction.objects.filter(public=True).update(public=False)

class Command(BaseCommand):
    def add_arguments(self, parser):
        parser.add_argument('filename', help='Path to csv file.')

    def handle(self, *args, **kwargs):
        fn = kwargs['filename']
        with open(fn) as f:
            reader = csv.Reader(f)
            reader.next() # skip headers
        
            for row in reader:
                ref = row[0]
        
                PerformedReaction.objects.filter(reference=ref).update(public=True)
        
