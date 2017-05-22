"""Recalculate the descriptors for all compounds and reactions."""
from django.core.management.base import BaseCommand
from DRP.models import Reaction, Compound
from django import db
from django.conf import settings
import logging
import importlib
from django.db import transaction

molDescriptorPlugins = [importlib.import_module(plugin) for
                        plugin in settings.MOL_DESCRIPTOR_PLUGINS]

rxnDescriptorPlugins = [importlib.import_module(plugin) for
                        plugin in settings.RXN_DESCRIPTOR_PLUGINS]

logger = logging.getLogger('DRP.management')


def calculate_descriptors(queryset, descriptorPlugins, verbose=False, plugins=None, **kwargs):
    """Helper function for descriptor calculation."""
    if verbose:
        logger.info(
            "Calculating descriptors for {} objects".format(queryset.count()))
    for plugin in descriptorPlugins:
        if plugins is None or plugin.__name__ in plugins:
            if verbose:
                logger.info("Calculating for plugin: {}".format(plugin))
            plugin.calculate_many(queryset, verbose=verbose, **kwargs)
            if verbose:
                logger.info("Done with plugin: {}\n".format(plugin))


class Command(BaseCommand):
    """Recalculate the descriptors for all compounds and reactions."""

    help = 'Recalculate the descriptors for all compounds and reactions.'

    def add_arguments(self, parser):
        """Add arguments for the parser."""
        parser.fromfile_prefix_chars = '@'
        parser.epilog = "Prefix arguments with '@' to specify a file containing newline-separated values for that argument. e.g.'-w @whitelist_headers.txt' to pass multiple descriptors from a file as whitelist."

        parser.add_argument('start', type=int, default=0, nargs='?',
                            help='pk of starting point. Indicates compound pk unless --reactions is specified')
        parser.add_argument('--count', type=int, default=5000, nargs='?',
                            help='Number of objects for this script to calculate for.')
        parser.add_argument('-e', '--error-level', nargs='?', default=0, const=3, type=int,
                            help='Make warnings errors instead. '
                                 '0 leaves python default settings '
                                 '(or whatever settings are specified by the command line flags when calling python). '
                                 '1 makes only RuntimeWarnings errors. '
                                 '2 makes Runtime and User Warnings errors. '
                                 '3 makes all warnings errors.')
        parser.add_argument('-p', '--plugins', nargs='+',
                            help='Plugins to use (default all).')
        parser.add_argument('-w', '--whitelist', nargs='+',
                            help='One or more descriptor headers to calculate from specified plugins (default all for given plugins).')
        group = parser.add_mutually_exclusive_group(required=False)
        group.add_argument('-r', '--reactions', '--rxns', action='store_true',
                           help='Calculate descriptors for reactions only.')
        group.add_argument('-c', '--compounds', action='store_true',
                           help='Calculate descriptors for compounds only.')
        group.add_argument('--include-invalid', action='store_true',
                           help='Calculate descriptors for invalid reactions also.')
        group.add_argument('--include-non-performed', action='store_true',
                           help='Calculate descriptors for non-performed reactions also.')
        group.add_argument('--only-dirty', action='store_true',
                           help='Calculate descriptors only for those objects which have been flagged for calculation.')

    def handle(self, *args, **kwargs):
        """Handle the function call."""
        verbose = (kwargs['verbosity'] > 0)
        only_reactions = kwargs['reactions']
        print(only_reactions)
        only_compounds = kwargs['compounds']
        start = kwargs['start']
        whitelist = kwargs['whitelist']
        plugins = kwargs['plugins']
        include_invalid = kwargs['include_invalid']
        include_non_performed = kwargs['include_non_performed']
        only_dirty = kwargs['only_dirty']
        print(only_dirty)
        limit = kwargs['count']

        if whitelist is not None:
            # just a little optimization
            whitelist = set(whitelist)

        if kwargs['error_level'] == 1:
            warnings.simplefilter('error', RuntimeWarning)
        if kwargs['error_level'] == 2:
            warnings.simplefilter('error', UserWarning)
            warnings.simplefilter('error', RuntimeWarning)
        if kwargs['error_level'] == 3:
            warnings.simplefilter('error')

        if not only_reactions:
            with transaction.atomic():
                compounds = Compound.objects.order_by('pk').filter(
                    pk__gte=start).exclude(calculating=True)
                logger.debug('Compounds count is {}'.format(compounds.count()))
                if only_dirty:
                    compounds = compounds.objects.filter(dirty=True)
                compounds = compounds[:limit]
                # This hits our database again, but we have to because slices
                # can't be updated and we need to call these specific reactions
                # back.'
                compounds = Compound.objects.filter(
                    id__in=(compound.id for compound in compounds))
                compounds.update(calculating=True)
            logger.debug('Compounds count is {}'.format(compounds.count()))
            while compounds.count() > 1:
                try:
                    calculate_descriptors(compounds, molDescriptorPlugins,
                                          verbose=verbose, whitelist=whitelist, plugins=plugins)
                except Exception as e:
                    compounds.update(calculating=False)
                    raise e
                with transaction.atomic():
                    compounds = compounds.all()  # Refresh the queryset
                    compounds.filter(recalculate=False).update(
                        dirty=False, calculating=False)
                    compounds = compounds.filter(recalculate=True)
                    compounds.update(recalculate=False)
        if not only_compounds:
            with transaction.atomic():
                reactions = Reaction.objects.order_by(
                    'pk').exclude(calculating=True)
                print(reactions.count())
#                reactions = reactions.exclude(compounds__dirty=True)
                print(reactions.count())
                if only_dirty:
                    reactions = reactions.objects.filter(dirty=True)
                if only_reactions:
                    reactions = reactions.filter(pk__gte=start)
                    print(reactions.count())
                if not include_invalid:
                    reactions = reactions.exclude(
                        performedreaction__valid=False)
                if not include_non_performed:
                    reactions = reactions.exclude(performedreaction=None)
                reactions = reactions[:limit]
                reactions = Reaction.objects.filter(
                    id__in=(reaction.id for reaction in reactions))
                reactions.update(calculating=True)
                print(reactions.count())
            while reactions.count() > 1:
                try:
                    calculate_descriptors(reactions, rxnDescriptorPlugins,
                                          verbose=verbose, whitelist=whitelist, plugins=plugins)
                except Exception as e:
                    reactions.update(calculating=False)
                    raise e
                with transaction.atomic():
                    reactions = reactions.all()  # refresh the qs
                    reactions.filter(recalculate=False).update(
                        dirty=False, calculating=False)
                    reactions = reactions.filter(recalculate=True)
                    reactions.update(recalculate=False)
