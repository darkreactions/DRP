"""Miniaturise the current database."""
from django.core.management.base import BaseCommand
from DRP.models import Reaction, Descriptor, Compound, CatRxnDescriptorValue, CatMolDescriptorValue, ModelContainer, DataSetRelation, DataSet
from django import db


class Command(BaseCommand):
    """Miniaturise the current database."""
    help = 'Miniaturize the current database.'

    def add_arguments(self, parser):
        """Add arguments to the parser."""
        parser.add_argument('-n', '--number', type=int, default=200,
                            help='Number of reactions to include.')
        parser.add_argument('-b', '--balanced', action='store_true',
                            help='Whether to have equal number of successful and failed reactions.')
        parser.add_argument('--valid-only', action='store_true',
                            help='Whether to only include valid reactions.')

    def handle(self, *args, **kwargs):
        """Handle the command call."""
        ModelContainer.objects.all().delete()
        DataSetRelation.objects.all().delete()
        DataSet.objects.all().delete()
        CatRxnDescriptorValue.objects.exclude(descriptor__calculatorSoftware='manual').delete()
        CatMolDescriptorValue.objects.exclude(descriptor__calculatorSoftware='manual').delete()
        Descriptor.objects.exclude(calculatorSoftware='manual').delete()

        rxns = Reaction.objects.all()

        if kwargs['valid_only']:
            rxns = rxns.filter(performedreaction__valid=True)
        if kwargs['balanced']:
            true_num = kwargs['number'] / 2
            false_num = kwargs['number'] - true_num

            true_rxns = rxns.filter(boolrxndescriptorvalue__descriptor__heading='boolean_crystallisation_outcome', boolrxndescriptorvalue__value=True).distinct().order_by('-pk')
            print true_rxns.count()
            cutoff_pk = true_rxns[true_num].pk
            true_rxns_to_delete = true_rxns.filter(pk__lte=cutoff_pk)
            print true_rxns_to_delete.count()
            assert(true_rxns.count() - true_rxns_to_delete.count() == true_num)
            true_rxns_to_delete.delete()

            false_rxns = rxns.filter(boolrxndescriptorvalue__descriptor__heading='boolean_crystallisation_outcome', boolrxndescriptorvalue__value=False).distinct().order_by('-pk')
            print false_rxns.count()
            cutoff_pk = false_rxns[false_num].pk
            false_rxns_to_delete = false_rxns.filter(pk__lte=cutoff_pk)
            print false_rxns_to_delete.count()
            assert(false_rxns.count() - false_rxns_to_delete.count() == false_num)
            false_rxns_to_delete.delete()

            rxns.exclude(boolrxndescriptorvalue__descriptor__heading='boolean_crystallisation_outcome').delete()
            rxns.filter(boolrxndescriptorvalue__descriptor__heading='boolean_crystallisation_outcome', boolrxndescriptorvalue__value=None).delete()
        else:
            rxns = rxns.order_by('-pk')
            cutoff_pk = rxns[kwargs['number']]
            rxns.filter(pk__lte=cutoff_pk).delete()

        print Reaction.objects.all().count()

        Compound.objects.filter(compoundquantity=None).delete()
