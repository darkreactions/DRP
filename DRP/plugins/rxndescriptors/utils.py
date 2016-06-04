'''A utilities module for helping with reaction descriptor plugins'''
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
            for k, v in descDict.items():
                args = v.copy()
                del args['type']
                args['heading'] = k
                fetchArgs = {k: v for k, v in args.items() if k in ('heading', 'calculatorSoftware', 'calculatorSoftwareVersion')}
                if v['type'] == 'num':
                    try:
                        self.internalDict[k] = DRP.models.NumRxnDescriptor.objects.filter(**fetchArgs).update(**args)
                        self.internalDict[k] = DRP.models.NumRxnDescriptor.objects.get(**fetchArgs)
                    except DRP.models.NumRxnDescriptor.DoesNotExist:
                        self.internalDict[k] = DRP.models.NumRxnDescriptor.objects.get_or_create(**args)[0]
                elif v['type'] == 'bool':
                    try:
                        self.internalDict[k] = DRP.models.BoolRxnDescriptor.objects.filter(**fetchArgs).update(**args)
                        self.internalDict[k] = DRP.models.BoolRxnDescriptor.objects.get(**fetchArgs)
                    except DRP.models.BoolRxnDescriptor.DoesNotExist:
                        self.internalDict[k] = DRP.models.BoolRxnDescriptor.objects.get_or_create(**args)[0]
                elif v['type'] == 'ord':
                    try:
                        self.internalDict[k] = DRP.models.OrdRxnDescriptor.objects.filter(**fetchArgs).update(**args)
                        self.internalDict[k] = DRP.models.OrdRxnDescriptor.objects.get(**fetchArgs)
                    except DRP.models.OrdRxnDescriptor.DoesNotExist:
                        self.internalDict[k] = DRP.models.OrdRxnDescriptor.objects.get_or_create(**args)[0]
                elif v['type'] == 'cat':
                    del args['permittedValues']
                    try:
                        DRP.models.CatRxnDescriptor.objects.filter(**fetchArgs).update(**args)
                        self.internalDict[k] = DRP.models.CatRxnDescriptor.objects.get(**fetchArgs)
                    except DRP.models.CatRxnDescriptor.DoesNotExist:
                        self.internalDict[k] = DRP.models.CatRxnDescriptor.objects.get_or_create(**args)[0]
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

    def __getattr__(self, name):
        """Deals with all names that are not defined explicitly by passing them to the internal dictionary (after initialising it).""""
        def _pass_attr(*args, **kwargs):
            self.initialise(self.descDict)
            return getattr(self.internalDict, name)(*args, **kwargs)

        return _pass_attr


def setup(descDict):
    return LazyDescDict(descDict)
