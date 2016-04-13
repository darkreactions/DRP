#!/usr/bin/env python
import django
django.setup()
from DRP.models import ModelContainer, NumRxnDescriptor
from sys import argv
import os
from django.db.models import Count

descriptor_directory = 'descs'

desc_files = ['small.dsc']
splitter = 'MutualInfoSplitter'
splitterOptions = '"num_splits": 15'
modelVisitorOptions = ['"BCR": true', '"BCR": false']
modelVisitorTools = ['J48', 'KNN', 'LogisticRegression', 'NaiveBayes', 'RandomForest', 'SVM_PUK']

containers = ModelContainer.objects.filter(built=True)

containers = containers.filter(splitterOptions__contains=splitterOptions)

for descriptor_file in desc_files:
    descriptor_file_path = os.path.join(descriptor_directory, descriptor_file)
    with open(descriptor_file_path) as f:
        headers = [l.strip() for l in f if l.strip()]

    numDescs = NumRxnDescriptor.objects.filter(heading__in=headers)
    num_descs = numDescs.count()

    if numDescs.count() != len(headers):
        raise RuntimeError("Did not find correct number of descriptors")

    #print numDescs
    conts = containers

    conts = conts.filter(numRxnDescriptors__in=numDescs).distinct()
    #print conts.count()

    conts = conts.annotate(num_descs=Count('numRxnDescriptors'))
    #for cont in conts:
        #print cont.num_descs, cont.numRxnDescriptors.all()
    conts = conts.filter(num_descs=num_descs)
    #print conts.count()
    
    for tool in modelVisitorTools:
        model_conts = conts.filter(modelVisitorTool=tool)
        print model_conts.count()
        for options in modelVisitorOptions:
            option_conts = model_conts.filter(modelVisitorOptions__contains=options)
            print descriptor_file, tool, options, option_conts.count()
