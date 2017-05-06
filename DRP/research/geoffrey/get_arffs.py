"""Get the arff files from a ModelContainer"""
import django
django.setup()
import os
from sys import argv
import shutil
from DRP.models import ModelContainer, StatsModel
from DRP.ml_models.model_visitors.weka import SVM_PUK
from itertools import chain

id = argv[1]
dest = argv[2]

c = ModelContainer.objects.get(id=id)

for i, sm in enumerate(c.statsmodel_set.all()):
    mv = SVM_PUK(statsModel=sm)
    fn = str(sm.inputFile)
    shutil.copyfile(fn, os.path.join(dest, str(i+1) + '_train_' + fn.split('/')[-1]))
    descriptorHeaders = [d.csvHeader for d in chain(c.descriptors, c.outcomeDescriptors)]

    assert(sm.testSets.all().count() == 1)
    testSet = sm.testSets.all()[0]
    test_arff = mv._prepareArff(testSet.reactions.all(), descriptorHeaders, verbose=True)
    shutil.move(test_arff, os.path.join(dest, str(i+1) + '_test_' + test_arff.split('/')[-1]))
