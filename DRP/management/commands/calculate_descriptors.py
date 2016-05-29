from django.core.management.base import BaseCommand
from DRP.models import Reaction, Compound
from django import db
import warnings

class Command(BaseCommand):
    help = 'Recalculate the descriptors for all compounds and reactions.'

    def add_arguments(self, parser):
        parser.fromfile_prefix_chars = '@'
        parser.epilog = "Prefix arguments with '@' to specify a file containing newline-separated values for that argument. e.g.'-w @whitelist_headers.txt' to pass multiple descriptors from a file as whitelist."

        parser.add_argument('start', type=int, default=0, nargs='?',
                    help='pk of starting point. Indicates compound pk unless --reactions is specified')
        parser.add_argument('-e', '--error-level', nargs='?', default=0, const=3, type=int,
                            help='Make warnings errors instead. '
                                 '0 leaves python default settings '
                                 '(or whatever settings are specified by the command line flags when calling python). '
                                 '1 makes only RuntimeWarnings errors. '
                                 '2 makes Runtime and User Warnings errors. '
                                 '3 makes all warnings errors.')
        parser.add_argument('-w', '--whitelist', nargs='+',
                            help='One or more descriptor headers to calculate')
        group = parser.add_mutually_exclusive_group(required=False)
        group.add_argument('-r', '--reactions', '--rxns', action='store_true',
                            help='Calculate descriptors for reactions only.')
        group.add_argument('-c', '--compounds', action='store_true',
                            help='Calculate descriptors for compounds only.')

    def handle(self, *args, **kwargs):
        verbose = (kwargs['verbosity'] > 0)
        only_reactions = kwargs['reactions']
        only_compounds = kwargs['compounds']
        start = kwargs['start']
        whitelist = kwargs['whitelist']

        if kwargs['error_level'] == 1:
            warnings.simplefilter('error', RuntimeWarning)
        if kwargs['error_level'] == 2:
            warnings.simplefilter('error', UserWarning)
            warnings.simplefilter('error', RuntimeWarning)
        if kwargs['error_level'] == 3:
            warnings.simplefilter('error')

        if not only_reactions:
            Compound.objects.order_by('pk').filter(pk__gte=start).calculate_descriptors(verbose=verbose, whitelist=whitelist)
        if not only_compounds:
            if only_reactions:
                Reaction.objects.order_by('pk').filter(pk__gte=start).calculate_descriptors(verbose=verbose, whitelist=whitelist)
            else:
                Reaction.objects.order_by('pk').calculate_descriptors(verbose=verbose, whitelist=whitelist)
