from django.core.management.base import BaseCommand
from DRP.models import Reaction, Compound
from django import db
import warnings

class Command(BaseCommand):
    help = 'Recalculate the descriptors for all compounds and reactions.'

    def add_arguments(self, parser):
        parser.add_argument('-e', '--error', action='store_true',
                            help='Make warnings errors instead')
        group = parser.add_mutually_exclusive_group(required=False)
        group.add_argument('-r', '--reactions', '--rxns', action='store_true',
                            help='Calculate descriptors for reactions only.')
        group.add_argument('-c', '--compounds', action='store_true',
                            help='Calculate descriptors for compounds only.')

    def handle(self, *args, **kwargs):
        verbose = (kwargs['verbosity'] > 0)
        only_reactions = kwargs['reactions']
        only_compounds = kwargs['compounds']
        if kwargs['error']:
            warnings.simplefilter('error')
        if not only_reactions:
            Compound.objects.all().calculate_descriptors(verbose=verbose)
        if not only_compounds:
            Reaction.objects.all().calculate_descriptors(verbose=verbose)
