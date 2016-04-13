#!/usr/bin/env python
import django
django.setup()
from DRP.models import ModelContainer, NumRxnDescriptor, BoolRxnDescriptor
from sys import argv
import os
from django.db.models import Count
import build_model

descriptor_directory = 'legacy_tests/final_descs/use/'


desc_files = ['CFS_new_legacy_noCA_nonZeroVariance.dsc']
splitter = 'MutualInfoSplitter'
splitterOptions = '"num_splits": 15'
modelVisitorOptions = ['"BCR": true', '"BCR": false']
modelVisitorTools = ['J48', 'KNN', 'LogisticRegression', 'NaiveBayes', 'RandomForest', 'SVM_PUK']

containers = ModelContainer.objects.filter(built=True).annotate(num_numDescs=Count('numRxnDescriptors', distinct=True)).annotate(num_boolDescs=Count('boolRxnDescriptors', distinct=True))

#containers = containers.filter(splitterOptions__contains=splitterOptions)

for descriptor_file in desc_files:
    descriptor_file_path = os.path.join(descriptor_directory, descriptor_file)
    with open(descriptor_file_path) as f:
        headers = [l.strip() for l in f if l.strip()]

    boolDescs = BoolRxnDescriptor.objects.filter(heading__in=headers)
    numDescs = NumRxnDescriptor.objects.filter(heading__in=headers)
    num_numDescs = numDescs.count()
    num_boolDescs = boolDescs.count()


    if numDescs.count() != len(headers):
        missing_descs = build_model.missing_descriptors(headers)
        raise RuntimeError("Did not find correct number of descriptors. Unable to find: {}".format(missing_descs))

    conts = containers
    if numDescs:
        conts = conts.filter(numRxnDescriptors__in=numDescs).distinct()
    if boolDescs:
        conts = conts.filter(boolRxnDescriptors__in=boolDescs).distinct()

    conts = conts.filter(num_numDescs=num_numDescs).filter(num_boolDescs=num_boolDescs)

    for tool in modelVisitorTools:
        model_conts = conts.filter(modelVisitorTool=tool)
        print model_conts.count()
        for options in modelVisitorOptions:
            option_conts = model_conts.filter(modelVisitorOptions__contains=options)
            print descriptor_file, tool, options, option_conts.count()
