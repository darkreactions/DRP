#!/usr/bin/env python
import django
from django.conf import settings
import importlib
from build_model import prepare_build_display_model

django.setup()
response_headers = ["boolean_crystallisation_outcome"]
predictor_headers = ["reaction_temperature"]
tools_to_test = (
                "SVM_PUK",
                "KNN",
                "NaiveBayes",
                "J48",
                "RandomForest",
                "LogisticRegression",
                )

visitorModules = {library:importlib.import_module(settings.STATS_MODEL_LIBS_DIR + "."+ library) for library in settings.STATS_MODEL_LIBS}

for modelVisitorLibrary, module in visitorModules.items():
    #for modelVisitorTool in module.tools:
    for modelVisitorTool in tools_to_test:
        # This only works for classifiers
        #if modelVisitorTool in tools_to_test:
        for splitter in settings.REACTION_DATASET_SPLITTERS:
            print modelVisitorLibrary, modelVisitorTool, splitter
            prepare_build_display_model(predictor_headers=predictor_headers, response_headers=response_headers, modelVisitorLibrary=modelVisitorLibrary, modelVisitorTool=modelVisitorTool, splitter=splitter, verbose=True)
    
    
