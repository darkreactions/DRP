
from django import forms
from django.forms.models import modelformset_factory
from DRP.models import PerformedReaction, OrdRxnDescriptorValue, CompoundQuantity
from DRP.models import NumRxnDescriptorValue, BoolRxnDescriptorValue, CatRxnDescriptorValue
from DRP.models import NumRxnDescriptor, BoolRxnDescriptor, CatRxnDescriptor, OrdRxnDescriptor

def descriptorValueFormFactory(modelClass, descriptorClass):

  class DescriptorValueForm(forms.ModelForm):

    class Meta:
      model=modelClass
      fields=('descriptor', 'value')

    def __init__(self, *args, **kwargs):
      super(DescriptorValueForm, self).__init__(*args, **kwargs)
      self.fields['descriptor'].queryset = descriptorClass.objects.filter(calculatorSoftware='manual')

  return DescriptorValueForm

factory = descriptorValueFormFactory

NumRxnDescValForm = factory(NumRxnDescriptorValue, NumRxnDescriptor)
OrdRxnDescValForm = factory(OrdRxnDescriptorValue, OrdRxnDescriptor)
BoolRxnDescValForm = factory(BoolRxnDescriptorValue, BoolRxnDescriptor)
CatRxnDescValForm = factory(CatRxnDescriptorValue, CatRxnDescriptor) 

NumRxnDescValFormSet = modelformset_factory(NumRxnDescriptorValue, form=NumRxnDescValForm, can_delete=True)
OrdRxnDescValFormSet = modelformset_factory(OrdRxnDescriptorValue, form=OrdRxnDescValForm, can_delete=True)
BoolRxnDescValFormSet = modelformset_factory(BoolRxnDescriptorValue, form=BoolRxnDescValForm, can_delete=True)
CatRxnDescValFormSet = modelformset_factory(CatRxnDescriptorValue, form=CatRxnDescValForm, can_delete=True)
