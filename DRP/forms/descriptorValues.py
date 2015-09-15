
from django import forms
from DRP.models import PerformedReaction, OrdRxnDescriptorValue, CompoundQuantity
from DRP.models import NumRxnDescriptorValue, BoolRxnDescriptorValue, CatRxnDescriptorValue
from DRP.models import NumRxnDescriptor, BoolRxnDescriptor, CatRxnDescriptor, OrdRxnDescriptor

def descriptorValueFormFactory(modelClass, descriptorClass):
'''A factory function for producing descriptor value form classes'''

  class DescriptorValueForm(forms.ModelForm):
    '''A form for a descriptor value'''

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
