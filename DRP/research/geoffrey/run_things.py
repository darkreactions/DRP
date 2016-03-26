import django
django.setup()
from DRP.models import PerformedReaction, DataSet, BoolRxnDescriptorValue

dsets = DataSet.objects.filter(name__contains='valid_legacy_rxns_nonzero_compound_239aed1e-e7c6-4215-9ba2-82a4287bb0f8')

sdsets = sorted([(int(d.name.split('_')[-1]), d) for d in dsets])
nsdsets = [(n,d) for n, d in sdsets if n%2==0]

true_weights = []
false_weights = []

for n, d in nsdsets:
    print d.name
    print n
    rxns = d.reactions.all()
    num_false = BoolRxnDescriptorValue.objects.filter(descriptor__heading='boolean_crystallisation_outcome', reaction__in=rxns, value=False).count()
    num_true = BoolRxnDescriptorValue.objects.filter(descriptor__heading='boolean_crystallisation_outcome', reaction__in=rxns, value=True).count()
    num_rxns = rxns.count()
    #print 'export TRAIN_NUM={}'.format(n)
    #print 'export TEST_NUM={}'.format(n+1)
    true_weights.append(num_rxns/float(num_true))
    false_weights.append(num_rxns/float(num_false))

print 'TRUE_WEIGHTS=({})'.format(' '.join(['"{}"'.format(w) for w in true_weights]))
print 'FALSE_WEIGHTS=({})'.format(' '.join(['"{}"'.format(w) for w in false_weights]))
