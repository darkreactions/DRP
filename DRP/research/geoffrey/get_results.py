#!/usr/bin/env python
import django
django.setup()
from DRP.models import ModelContainer, NumRxnDescriptor, BoolRxnDescriptor
from sys import argv
import os
from django.db.models import Count
import build_model

descriptor_directory = 'legacy_tests/final_descs/use/'


desc_files = ['new_leak_slowcool_group_period_valence_nonZeroVariance.dsc']
splitter = 'MutualInfoSplitter'
splitterOptions = '"num_splits": 15'
modelVisitorOptions = ['"BCR": true', '"BCR": false']
modelVisitorTools = ['J48', 'KNN', 'LogisticRegression', 'NaiveBayes', 'RandomForest', 'SVM_PUK']

containers = ModelContainer.objects.filter(built=True).annotate(num_numDescs=Count('numRxnDescriptors', distinct=True)).annotate(num_boolDescs=Count('boolRxnDescriptors', distinct=True))

print "{} built model containers".format(containers.count())
containers = containers.filter(splitter=splitter).filter(splitterOptions__contains=splitterOptions)

print "{} with correct splitter".format(containers.count())

for descriptor_file in desc_files:
    print "Using {}".format(descriptor_file)
    descriptor_file_path = os.path.join(descriptor_directory, descriptor_file)
    with open(descriptor_file_path) as f:
        headers = [l.strip() for l in f if l.strip()]

    boolDescs = BoolRxnDescriptor.objects.filter(heading__in=headers)
    numDescs = NumRxnDescriptor.objects.filter(heading__in=headers)
    num_numDescs = numDescs.count()
    num_boolDescs = boolDescs.count()

    if num_numDescs + num_boolDescs != len(headers):
        missing_descs = []
        for header in headers:
            if not NumRxnDescriptor.objects.filter(heading=header).exists() and not BoolRxnDescriptor.objects.filter(heading=header).exists():
                missing_descs.append(header)
        raise RuntimeError("Did not find correct number of descriptors. Unable to find:{}".format(missing_descs))

    conts = containers
    conts = conts.filter(num_numDescs=num_numDescs).filter(num_boolDescs=num_boolDescs)
    print "{} containers with appropriate number of descriptors".format(conts.count())
    
    for tool in modelVisitorTools:
        model_conts = conts.filter(modelVisitorTool=tool)
        print model_conts.count()
        for options in modelVisitorOptions:
            option_conts = model_conts.filter(modelVisitorOptions__contains=options)
            print descriptor_file, tool, options, option_conts.count()
            for cont in option_conts:
                print cont.id
