#!/usr/bin/env python
import django
django.setup()
from DRP.models import ModelContainer, NumRxnDescriptor, BoolRxnDescriptor
from sys import argv
import os
from django.db.models import Count
from DRP import utils
import json
import csv
from math import sqrt

descriptor_directory = 'legacy_tests/final_descs/use/'

# dictionary from file name to descriptor set specification
# 
desc_files = {
                'CFS_legacy_mw_noCA_nonZeroVariance.dsc' : {'Base Features':'legacy', 'ChemAxon':False, 'Feature Selection':'CFS'},
                'CFS_legacy_mw_noPSA_nonZeroVariance.dsc' : {'Base Features':'legacy', 'ChemAxon':True, 'Feature Selection':'CFS'},
                'CFS_new_leak_slowcool_group_period_valence_nonZeroVariance.dsc' : {'Base Features':'new', 'ChemAxon':True, 'Feature Selection':'CFS'},
                'CFS_new_legacy_bothCA_noPSA_nonZeroVariance.dsc' : {'Base Features':'new + legacy', 'ChemAxon':'Both', 'Feature Selection':'CFS'},
                'CFS_new_legacy_legCA_noPSA_nonZeroVariance.dsc' : {'Base Features':'new + legacy', 'ChemAxon':'Legacy', 'Feature Selection':'CFS'},
                'CFS_new_legacy_newCA_nonZeroVariance.dsc' : {'Base Features':'new + legacy', 'ChemAxon':'New', 'Feature Selection':'CFS'},
                'CFS_new_legacy_noCA_nonZeroVariance.dsc' : {'Base Features':'new + legacy', 'ChemAxon':False, 'Feature Selection':'CFS'},
                'CFS_new_noCA_leak_slowcool_group_period_valence_nonZeroVariance.dsc' : {'Base Features':'new', 'ChemAxon':False, 'Feature Selection':'CFS'},
                'legacy_mw_noCA_nonZeroInfo.dsc' : {'Base Features':'legacy', 'ChemAxon':False, 'Feature Selection':'Postive Info'},
                'legacy_mw_noCA_nonZeroVariance.dsc' : {'Base Features':'legacy', 'ChemAxon':False, 'Feature Selection':'Postive Variance'},
                'legacy_mw_noPSA_nonZeroInfo.dsc' : {'Base Features':'legacy', 'ChemAxon':True, 'Feature Selection':'Postive Info'},
                'legacy_mw_noPSA_nonZeroVariance.dsc' : {'Base Features':'legacy', 'ChemAxon':True, 'Feature Selection':'Postive Variance'},
                'new_leak_slowcool_group_period_valence_nonZeroInfo.dsc' : {'Base Features':'new', 'ChemAxon':True, 'Feature Selection':'Postive Info'},
                'new_leak_slowcool_group_period_valence_nonZeroVariance.dsc' : {'Base Features':'new', 'ChemAxon':True, 'Feature Selection':'Postive Variance'},
                'new_legacy_bothCA_noPSA_nonZeroInfo.dsc' : {'Base Features':'new + legacy', 'ChemAxon':'Both', 'Feature Selection':'Postive Info'},
                'new_legacy_bothCA_noPSA_nonZeroVariance.dsc' : {'Base Features':'new + legacy', 'ChemAxon':'Both', 'Feature Selection':'Postive Variance'},
                'new_legacy_legCA_noPSA_nonZeroInfo.dsc' : {'Base Features':'new + legacy', 'ChemAxon':'Legacy', 'Feature Selection':'Postive Info'},
                'new_legacy_legCA_noPSA_nonZeroVariance.dsc' : {'Base Features':'new + legacy', 'ChemAxon':'Legacy', 'Feature Selection':'Postive Variance'},
                'new_legacy_newCA_nonZeroInfo.dsc' : {'Base Features':'new + legacy', 'ChemAxon':'New', 'Feature Selection':'Postive Info'},
                'new_legacy_newCA_nonZeroVariance.dsc' : {'Base Features':'new + legacy', 'ChemAxon':'New', 'Feature Selection':'Postive Variance'},
                'new_legacy_noCA_nonZeroInfo.dsc' : {'Base Features':'new + legacy', 'ChemAxon':False, 'Feature Selection':'Postive Info'},
                'new_legacy_noCA_nonZeroVariance.dsc' : {'Base Features':'new + legacy', 'ChemAxon':False, 'Feature Selection':'Postive Variance'},
                'new_noCA_leak_slowcool_group_period_valence_nonZeroInfo.dsc' : {'Base Features':'new', 'ChemAxon':False, 'Feature Selection':'Postive Info'},
                'new_noCA_leak_slowcool_group_period_valence_nonZeroVariance.dsc' : {'Base Features':'new', 'ChemAxon':False, 'Feature Selection':'Postive Variance'},
             }
splitter = 'RandomSplitter'
splitterOptions = json.dumps({"num_splits": 15})
BCR_options = [True, False]
modelVisitorTools = ['J48', 'KNN', 'LogisticRegression', 'NaiveBayes', 'RandomForest', 'SVM_PUK']

#with open('test.csv', 'w') as f:
    #headers = ['Base Features', 'ChemAxon', 'Feature Selection']
    #writer = csv.DictWriter(f, headers)
    #writer.writeheader()
    #for desc_file_spec in desc_files.values():
        #writer.writerow(desc_file_spec)
#exit()

def get_rows():
    rows = []
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
                row['Model'] = tool
                row['BCR Weighted'] = bcr_option
                
                if len(option_conts) != 1:
                    print "Was unable to find a unique model container matching given specification {}. Found {}".format(row, len(option_conts))
                    if len(option_conts) > 1:
                        print "Using container with largest pk"
                        option_conts.sort(key=lambda x: x.pk, reverse=True)
                    else:
                        print "Skipping this setup"
                        continue
                
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
    
                TP = true_guesses[True]
                FN = true_guesses[False]
                TN = false_guesses[False]
                FP = false_guesses[True]
                row['TP'] = TP
                row['FN'] = FN
                row['TN'] = TN
                row['FP'] = FP
                
                PP = TP+FP
                AP = TP+FN
                AN = TN+FP
                PN = TN+FN
                
                row['Accuracy'] = TP + TN
                row['BCR']  = (TP/AP + TN/AN)/2
                row['Matthews'] = (TP*TN - FP*FN)/sqrt(PP*AP*AN*PN) if 0 not in [PP, AP, AN, PN] else 0.0
                rows.append(row)
    return rows
    
if __name__ == '__main__':
    rows = get_rows()

    headers = ['Base Features', 'ChemAxon', 'Feature Selection', 'BCR Weighted', 'Model', 'TP', 'FP', 'FN', 'TN', 'Accuracy', 'BCR', 'Matthews']
    csv_fn = argv[1]

    with open(csv_fn, 'w') as f:
        writer = csv.DictWriter(f, headers)
        writer.writeheader()
        writer.writerows(rows)
