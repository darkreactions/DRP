#!/usr/bin/env python

from DRP.models import PerformedReaction, ModelContainer, Descriptor, rxnDescriptorValues
from DRP.models.rxnDescriptorValues import OrdRxnDescriptorValue, NumRxnDescriptorValue, BoolRxnDescriptorValue, CatRxnDescriptorValue
from django.db.models import Q
import operator

def build_model():
  reactions = PerformedReaction.objects.all()
  
  container = ModelContainer("weka", "SVM_PUK_basic", splitter="KFoldSplitter",
                             reactions=reactions)
  container.save()

  headers = ["reaction_temperature"]

  predictors = get_descriptors_by_header(headers)
  responses = Descriptor.objects.filter(heading="boolean_crystallisation_outcome")

  container.build(predictors, responses)

  print container.summarize()

def get_descriptors_by_header(headers):
  Qs = [Q(heading=header) for header in headers]
  return Descriptor.objects.filter(reduce(operator.or_, Qs))

if __name__=='__main__':
  build_model()
