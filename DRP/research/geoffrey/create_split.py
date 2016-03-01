#!/usr/bin/env python
import django
from DRP.models import PerformedReaction, MetricContainer, Descriptor, rxnDescriptorValues, DataSet
from DRP.ml_models.splitters.SingleSplitter import Splitter as SingleSplitter

def build_training_test_set(name, response_headers=["boolean_crystallisation_outcome"]):
    reactions = PerformedReaction.objects.filter(valid=True)
    reactions = reactions.exclude(ordrxndescriptorvalue__in=rxnDescriptorValues.OrdRxnDescriptorValue.objects.filter(descriptor__heading__in=response_headers, value=None))
    reactions = reactions.exclude(boolrxndescriptorvalue__in=rxnDescriptorValues.BoolRxnDescriptorValue.objects.filter(descriptor__heading__in=response_headers, value=None))
    reactions = reactions.exclude(catrxndescriptorvalue__in=rxnDescriptorValues.CatRxnDescriptorValue.objects.filter(descriptor__heading__in=response_headers, value=None))

    splitter = SingleSplitter(name)

    splits = splitter.split(reactions, verbose=True)



if __name__ == '__main__':
    django.setup()
    build_training_test_set("hello_testing_datasets_again")
