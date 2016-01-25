#!/usr/bin/env python

from DRP.models import PerformedReaction, ModelContainer, Descriptor, rxnDescriptorValues
from DRP.models.rxnDescriptorValues import OrdRxnDescriptorValue, NumRxnDescriptorValue, BoolRxnDescriptorValue, CatRxnDescriptorValue

def build_model():
  reactions = PerformedReaction.objects.all()
  
  container = ModelContainer("weka", "SVM", splitter="KFoldSplitter",
                             reactions=reactions)
  container.save()

  predictors = Descriptor.objects.filter(heading="reaction_temperature")[:100]
  responses = Descriptor.objects.filter(heading="boolean_crystallisation_outcome")[:100]

  container.build(predictors, responses)

  print container.summarize()

if __name__=='__main__':
  build_model()
