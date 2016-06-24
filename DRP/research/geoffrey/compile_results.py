import django
django.setup()

from DRP.models import ModelContainer, Descriptor, NumRxnDescriptor, BoolRxnDescriptor, OrdRxnDescriptor, CatRxnDescriptor
import os.path

descriptor_directory = 'descs'
descriptor_filenames = ['base_noCA.dsc',
                        ]


model_visitors = ['SVM_PUK',
                  'J48',
                  'KNN',
                  'LogisticRegression',
                  'NaiveBayes',
                  'RandomForest',
                  ]

splitters = ['RandomSplitter',
             'ExploratorySplitter',
             ]

response = BoolRxnDescriptor.objects.filter(heading='boolean_crystallisation_outcome')

built_containers = ModelContainer.objects.filter(built=True)

for fn in descriptor_filenames:
    path = os.path.join(descriptor_directory, fn)
    with open(path) as f:
        headings = [h.strip() for h in f.readlines() if h]

    num_predictors = NumRxnDescriptor.objects.filter(heading__in=headings)
    bool_predictors = BoolRxnDescriptor.objects.filter(heading__in=headings)
    for visitor in model_visitors:
        containers = built_containers.filter(outcomeBoolRxnDescriptors=response, numRxnDescriptors=num_predictors, boolRxnDescriptors=bool_predictors, modelVisitorTool=visitor).distinct()
        for c in containers:
            print c
