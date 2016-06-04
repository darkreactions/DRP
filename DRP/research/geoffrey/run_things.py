#!/usr/bin/env python
import django
django.setup()
from DRP.models import ModelContainer, NumRxnDescriptor
from sys import argv
import os
from django.db.models import Count

containers = ModelContainer.objects.filter(built=True).only('numRxnDescriptors', 'id')
print containers.query
#exit()
print containers.count()
descriptor_directory = 'descs'

desc_files = ['small.dsc']
modelVisitorOptions = ['{"BCR": true}', '{"BCR": false}']
modelVisitorTools = ['J48', 'KNN', 'LogisticRegression', 'NaiveBayes', 'RandomForest', 'SVM_PUK']


for descriptor_file in desc_files:
    descriptor_file_path = os.path.join(descriptor_directory, descriptor_file)
    with open(descriptor_file_path) as f:
        headers = [l.strip() for l in f if l.strip()]

    numDescs = list(NumRxnDescriptor.objects.filter(heading__in=headers))
    num_descs = len(numDescs) #numDescs.count()

    if num_descs != len(headers):
        raise RuntimeError("Did not find correct number of descriptors")

    #print numDescs
    conts = containers

    print numDescs
    #for c in conts:
        #print c.numRxnDescriptors.all()
    conts = conts.filter(numRxnDescriptors__in=numDescs)
    print conts.count()
    #for c in conts:
        #print c.num_descs
    #conts = conts
    conts = conts.annotate(num_descs=Count('numRxnDescriptors')).filter(num_descs=num_descs)
    print conts.count()
    print conts.query
