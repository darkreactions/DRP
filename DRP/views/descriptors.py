# this is the view file for the descriptors in add reaction

# import all of the things
from django.shortcuts import render
from django.views.generic import ListView
from DRP.models import BoolRxnDescriptor as BoolRxnDescriptorModel, CatRxnDescriptor as CatRxnDescriptorModel, NumRxnDescriptor as NumRxnDescriptorModel, OrdRxnDescriptor as OrdRxnDescriptorModel

# fetch descriptor details (from python, export to .json)


class BoolRxnDescriptor(ListView):

    def __init__(self, *args, **kwargs):
        self.queryset = BoolRxnDescriptorModel.objects.all()
        super(ListView, self).__init__(*args, **kwargs)

    def dispatch(self, request):
        return render(request, 'bool_rxn_descriptors.json', {'descriptors': self.queryset})


class NumRxnDescriptor(ListView):

    def __init__(self):
        self.des = NumRxnDescriptorModel.objects.all()

    def dispatch(self, request):
        return render(request, 'num_rxn_descriptors.json', {'descriptors': self.queryset})


class OrdRxnDescriptor(ListView):

    def __init__(self):
        self.des = OrdRxnDescriptorModel.objects.all()

    def dispatch(self, request):
        return render(request, 'ord_rxn_descriptors.json', {'descriptors': self.queryset})


class CatRxnDescriptor(ListView):

    def __init__(self):
        self.des = CatRxnDescriptorModel.objects.all()

    def dispatch(self, request):
        return render(request, 'cat_rxn_descriptors.json', {'descriptors': self.queryset})


# these find the descriptors and return them as a json object
#		as a list...
# tuples of descriptor and type

# not a .json file... .json object...

# look at the reactions.py
