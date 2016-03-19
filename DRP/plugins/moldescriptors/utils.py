'''A utilities module for helping with molecular descriptor plugins'''
import DRP
from django.db import transaction

class LazyDescDict(object):

    def __init__(self, descDict):
        self.internalDict = {}
        self.descDict = descDict
        self.initialised = False


    @transaction.atomic
    def initialise(self, descDict):
        if not self.initialised:
            for k,v in descDict.items():
                args = v.copy()
                del args['type']
                args['heading']=k
                fetchArgs = {k:v for k, v in args.items() if k in ('heading', 'calculatorSoftware', 'calculatorSoftwareVersion')}
                if v['type'] == 'num':
                    try:
                        self.internalDict[k] = DRP.models.NumMolDescriptor.objects.filter(**fetchArgs).update(**args)
                        self.internalDict[k] = DRP.models.NumMolDescriptor.objects.get(**fetchArgs)
                    except DRP.models.NumMolDescriptor.DoesNotExist:
                        self.internalDict[k] = DRP.models.NumMolDescriptor.objects.create(**args)
                elif v['type'] == 'bool':
                    try:
                        self.internalDict[k] = DRP.models.BoolMolDescriptor.objects.filter(**fetchArgs).update(**args)
                        self.internalDict[k] = DRP.models.BoolMolDescriptor.objects.get(**fetchArgs)
                    except DRP.models.BoolMolDescriptor.DoesNotExist:
                        self.internalDict[k] = DRP.models.BoolMolDescriptor.objects.create(**args)
                elif v['type'] == 'ord':
                    try:
                        self.internalDict[k] = DRP.models.OrdMolDescriptor.objects.filter(**fetchArgs).update(**args)
                        self.internalDict[k] = DRP.models.OrdMolDescriptor.objects.get(**fetchArgs)
                    except DRP.models.OrdMolDescriptor.DoesNotExist:
                        self.internalDict[k] = DRP.models.OrdMolDescriptor.objects.create(**args)
                elif v['type'] == 'cat':
                    del args['permittedValues']
                    try:
                        DRP.models.CatMolDescriptor.objects.filter(**fetchArgs).update(**args)
                        self.internalDict[k] = DRP.models.CatMolDescriptor.objects.get(**fetchArgs)
                    except DRP.models.CatMolDescriptor.DoesNotExist:
                        self.internalDict[k] = DRP.models.CatMolDescriptor.objects.create(**args)
                    for permittedValue in v['permittedValues']:
                        perm = DRP.models.CategoricalDescriptorPermittedValue.objects.get_or_create(value=permittedValue, descriptor=self.internalDict[k])[0]
                else:
                    raise RuntimeError("Invalid descriptor type provided")
        self.initialised = True

    def __len__(self):
        return len(self.internalDict)

    def __iter__(self):
        self.initialise(self.descDict)
        return iter(self.internalDict)

    def __getitem__(self, key):
        self.initialise(self.descDict)
        return self.internalDict[key]

    def __contains__(self, item):
        return item in self.internalDict


def setup(descDict):
    return LazyDescDict(descDict)
