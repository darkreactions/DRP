#!/usr/bin/env python
import django
django.setup()
from DRP.models import ModelContainer, NumRxnDescriptor, BoolRxnDescriptor
from sys import argv
import os
from django.db.models import Count
from DRP import utils
import json
import warnings

descriptor_directory = 'descs' #'legacy_tests/final_descs/use/'

# dictionary from file name to descriptor set specification
# 
desc_files = {'small_with_bool.dsc': {'base':'small', 'ChemAxon':False, 'Feature Selection':'Sample'}}
splitter = 'MutualInfoSplitter'
splitterOptions = json.dumps({"num_splits": 15})
BCR_options = [True, False]
modelVisitorTools = ['J48', 'KNN', 'LogisticRegression', 'NaiveBayes', 'RandomForest', 'SVM_PUK']

containers = ModelContainer.objects.filter(built=True).annotate(num_numDescs=Count('numRxnDescriptors', distinct=True)).annotate(num_boolDescs=Count('boolRxnDescriptors', distinct=True))

print "{} built model containers".format(containers.count())
containers = containers.filter(splitter=splitter).filter(splitterOptions=splitterOptions)

print "{} with correct splitter".format(containers.count())

for descriptor_file, desc_set_spec in desc_files.items():
    print "Using {}".format(descriptor_file)
    descriptor_file_path = os.path.join(descriptor_directory, descriptor_file)
    with open(descriptor_file_path) as f:
        headers = [l.strip() for l in f if l.strip()]

    boolDescs = BoolRxnDescriptor.objects.filter(heading__in=headers)
    numDescs = NumRxnDescriptor.objects.filter(heading__in=headers)
    num_numDescs = numDescs.count()
    num_boolDescs = boolDescs.count()

    setBoolDescs = set(boolDescs)
    setNumDescs = set(numDescs)

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
        for bcr_option in BCR_options:
            options = {'BCR': bcr_option}
            options_string = json.dumps(options)
            option_conts = model_conts.filter(modelVisitorOptions=options_string)
            option_conts = [c for c in option_conts if (set(c.numRxnDescriptors.all()) == setNumDescs and set(c.boolRxnDescriptors.all()) == setBoolDescs)]
            row = desc_set_spec.copy()
            row.update(options)
            row['model'] = tool
            if len(option_conts) != 1:
                print "Was unable to find a unique model container matching given specification {}. Found {}".format(row, len(option_conts))
            else:
                print "Found unique container"
                cont = option_conts[0]
                conf_tuples_lol = cont.getComponentConfusionMatrices()
                # we only care about the confusion matrix of the first descriptor
                confs = [conf_tuple_list[0][1] for conf_tuple_list in conf_tuples_lol]
                average_conf = utils.average_normalized_conf(confs)
                # And here we go into the part that only works for boolean descriptors
                # That's what I need right now and we're probably going to overhaul the
                # confusion matrix stuff soon anyway
                true_guesses = average_conf[True]
                false_guesses = average_conf[False]
    
                row['TP'] = true_guesses[True]
                row['FN'] = true_guesses[False]
                row['TN'] = false_guesses[True]
                row['FP'] = false_guesses[False]
                print row
    
