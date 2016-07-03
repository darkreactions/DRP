#!/usr/bin/env python
import django
django.setup()
from DRP.models import ModelContainer, NumRxnDescriptor, BoolRxnDescriptor
import sys
import os
from django.db.models import Count
from DRP import utils
import json
import csv
from math import sqrt

descriptor_directory = 'thesis_paper/descs'

# dictionary from file name to descriptor set specification
#
_desc_files = {
    # 'base_noCA.dsc': {'Descriptors': 'base', 'ChemAxon': False, 'Feature Selection': 'Variance+'},
    # 'base_CA.dsc': {'Descriptors': 'base', 'ChemAxon': True, 'Feature Selection': 'Variance+'},
    # 'orthogonal_noCA.dsc': {'Descriptors': 'orthogonal', 'ChemAxon': False, 'Feature Selection': 'Variance+'},
    # 'orthogonal_CA.dsc': {'Descriptors': 'orthogonal', 'ChemAxon': True, 'Feature Selection': 'Variance+'},
    # 'orthogonal_CA_noInorgCA.dsc': {'Descriptors': 'orthogonal', 'ChemAxon': 'Org only', 'Feature Selection': 'Variance+'},
    'orthogonal_plus_noCA.dsc': {'Descriptors': 'orthogonal+', 'ChemAxon': False, 'Feature Selection': 'Variance+'},
    # 'orthogonal_plus_CA.dsc': {'Descriptors': 'orthogonal+', 'ChemAxon': True, 'Feature Selection': 'Variance+'},
    # 'orthogonal_plus_CA_noInorgCA.dsc': {'Descriptors': 'orthogonal+', 'ChemAxon': 'Org only', 'Feature Selection': 'Variance+'},
}

desc_files = {}

for fn, attr in _desc_files.items():
    info_fn = fn[:-4] + '_info.dsc'
    CFS_fn = 'CFS_' + fn

    new_attr = attr.copy()
    new_attr['Feature Selection'] = 'Info+'
    desc_files[info_fn] = new_attr

    new_attr = attr.copy()
    new_attr['Feature Selection'] = 'CFS'
    desc_files[CFS_fn] = new_attr

desc_files.update(_desc_files)

splitter = 'ExploratorySplitter'
splitterOptions = json.dumps({"num_splits": 15})
BCR_options = [True, False]
modelVisitorTools = ['J48', 'KNN', 'LogisticRegression',
                     'NaiveBayes', 'RandomForest', 'SVM_PUK']

# with open('test.csv', 'w') as f:
#headers = ['Base Features', 'ChemAxon', 'Feature Selection']
#writer = csv.DictWriter(f, headers)
# writer.writeheader()
# for desc_file_spec in desc_files.values():
# writer.writerow(desc_file_spec)
# exit()


def get_rows():
    rows = []
    containers = ModelContainer.objects.filter(built=True).annotate(num_numDescs=Count(
        'numRxnDescriptors', distinct=True)).annotate(num_boolDescs=Count('boolRxnDescriptors', distinct=True))

    print "{} built model containers".format(containers.count())
    containers = containers.filter(splitter=splitter).filter(
        splitterOptions=splitterOptions)

    print "{} with correct splitter".format(containers.count())

    for descriptor_file, desc_set_spec in desc_files.items():
        print "Using {}".format(descriptor_file)
        descriptor_file_path = os.path.join(
            descriptor_directory, descriptor_file)
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
            raise RuntimeError(
                "Did not find correct number of descriptors. Unable to find:{}".format(missing_descs))

        conts = containers
        conts = conts.filter(num_numDescs=num_numDescs).filter(
            num_boolDescs=num_boolDescs)
        print "{} containers with appropriate number of descriptors".format(conts.count())

        for tool in modelVisitorTools:
            model_conts = conts.filter(modelVisitorTool=tool)
            for bcr_option in BCR_options:
                options = {'BCR': bcr_option}
                options_string = json.dumps(options)
                option_conts = model_conts.filter(
                    modelVisitorOptions=options_string)
                option_conts = [c for c in option_conts if (set(c.numRxnDescriptors.all(
                )) == setNumDescs and set(c.boolRxnDescriptors.all()) == setBoolDescs)]

                row = desc_set_spec.copy()
                row['Model'] = tool
                row['BCR Weighted'] = bcr_option

                if len(option_conts) != 1:
                    sys.stderr.write("Was unable to find a unique model container matching given specification {}. Found {}\n".format(row, len(option_conts)))
                    if len(option_conts) > 1:
                        sys.stderr.write("Using container with largest pk\n")
                        option_conts.sort(key=lambda x: x.pk, reverse=True)
                    else:
                        print sys.stderr.write("Skipping this setup\n")
                        continue

                cont = option_conts[0]
                print "{}\t{}".format(row, cont.id)
                conf_tuples_lol = cont.getComponentConfusionMatrices()
                # we only care about the confusion matrix of the first
                # descriptor
                confs = [conf_tuple_list[0][1]
                         for conf_tuple_list in conf_tuples_lol]
                average_conf = utils.average_normalized_conf(confs)

                # And here we go into the part that only works for boolean descriptors
                # That's what I need right now and we're probably going to overhaul the
                # confusion matrix stuff soon anyway
                true_guesses = average_conf[True]
                false_guesses = average_conf[False]

                TP = true_guesses[True]
                FN = true_guesses[False]
                TN = false_guesses[False]
                FP = false_guesses[True]
                row['TP'] = TP
                row['FN'] = FN
                row['TN'] = TN
                row['FP'] = FP

                PP = TP + FP
                AP = TP + FN
                AN = TN + FP
                PN = TN + FN

                row['Accuracy'] = TP + TN
                row['BCR'] = (TP / AP + TN / AN) / 2
                row['Matthews'] = (TP * TN - FP * FN) / sqrt(PP *
                                                             AP * AN * PN) if 0 not in [PP, AP, AN, PN] else 0.0
                rows.append(row)
    return rows

if __name__ == '__main__':
    rows = get_rows()

    headers = ['Descriptors', 'ChemAxon', 'Feature Selection', 'BCR Weighted',
               'Model', 'TP', 'FP', 'FN', 'TN', 'Accuracy', 'BCR', 'Matthews']
    csv_fn = sys.argv[1]

    with open(csv_fn, 'w') as f:
        writer = csv.DictWriter(f, headers)
        writer.writeheader()
        writer.writerows(rows)
