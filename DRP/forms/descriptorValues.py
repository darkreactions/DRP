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

class CatRxnDescValForm(forms.ModelForm):
    '''Special case because of the need to group permitted values'''
    
    value = forms.ChoiceField()

    class Meta:
        model=CatRxnDescriptorValue
        fields=('descriptor', 'value')

    def __init__(self, *args, **kwargs):
        super(CatRxnDescValForm, self).__init__(*args, **kwargs)
        descriptors = CatRxnDescriptor.objects.filter(calculatorSoftware='manual').prefetch_related('permittedValues')
        self.fields['descriptor'].queryset = descriptors
        self.fields['value'].choices = tuple(((descriptor.name, tuple((value.id, value.value) for value in descriptor.permittedValues.all())) for descriptor in descriptors))
