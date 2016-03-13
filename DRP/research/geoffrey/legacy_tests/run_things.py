import django
django.setup()
from DRP.models import PerformedReaction, DataSet
from django.db.models import Count

reactions = DataSet.objects.get(name='legacy_reactions').reactions.all()
#reactions = PerformedReaction.objects.all()
#reactions = reactions.filter(valid=True)

reactions = reactions.annotate(compound_count=Count('compounds'))
with open('legacy_reactions_allDescs.csv', 'wb') as f:
    reactions.toCsv(f, expanded=True)

#compound_numbers = [reactions.filter(compound_count=i).count() for i in xrange(10)]
#print compound_numbers


#rxn0compound = reactions.filter(compound_count=0)
#for rxn in rxn0compound:
    #print rxn
