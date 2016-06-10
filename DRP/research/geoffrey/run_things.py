#!/usr/bin/env python
import django
django.setup()
import DRP
from DRP.models import Descriptor
from DRP.models.querysets import MultiQuerySet


def rxn_descriptors():
    return MultiQuerySet(
                        NumRxnDescriptor.objects.all(),
                        BoolRxnDescriptor.objects.all(),
                        CatRxnDescriptor.objects.all(),
                        OrdRxnDescriptor.objects.all(),
                        )


if __name__ == '__main__':
    descs = rxn_descriptors()
    exclude_substring = ["_prediction_", "outcome", "rxnSpaceHash", "examplepy"]
    exclude_prefix = ["_", "transform"]

