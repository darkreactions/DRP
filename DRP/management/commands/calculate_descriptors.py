from django.core.management.base import BaseCommand
from DRP.models import Reaction, Compound
from django import db


class Command(BaseCommand):
    help = 'Recalculate the descriptors for all compounds and reactions.'

    def handle(self, *args, **kwargs):
        verbose = (kwargs['verbosity'] > 0)
        Compound.objects.all().calculate_descriptors(verbose=verbose)
        Reaction.objects.all().calculate_descriptors(verbose=verbose)
