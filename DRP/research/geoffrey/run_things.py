#!/usr/bin/env python
import django
django.setup()
from DRP.models import CompoundQuantity

for cq in CompoundQuantity.objects.filter(amount=0):
    print cq
    rxn = cq.reaction.performedreaction
    print rxn
    rxn.valid = False
    rxn.notes += 'Compound {} had amount {}. Amount changed None and reaction invalidated.'.format(cq.compound, cq.amount)
    rxn.save(calcDescriptors=False, invalidate_models=False)
    cq.amount = None
    cq.save(calcDescriptors=False, invalidate_models=False)
