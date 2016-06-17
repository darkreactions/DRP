"""Factories for descriptor value Forms."""
from django import forms
from DRP.models import PerformedReaction, OrdRxnDescriptorValue, CompoundQuantity, Reaction
from DRP.models import NumRxnDescriptorValue, BoolRxnDescriptorValue, CatRxnDescriptorValue
from DRP.models import NumRxnDescriptor, BoolRxnDescriptor, CatRxnDescriptor, OrdRxnDescriptor


def descriptorValueFormFactoryFactory(modelClass, descriptorClass):
    """A factory function for producing descriptor value form classes."""
    def descriptorValueFormFactory(reactionId):
        """A factory for descriptor value forms."""
        class DescriptorValueForm(forms.ModelForm):

            """A form for a descriptor value."""
    
            class Meta:
                model = modelClass
                fields = ('descriptor', 'value', 'reaction')
    
            def __init__(self, *args, **kwargs):
                """Lock the reaction to the specified value."""
                super(DescriptorValueForm, self).__init__(*args, **kwargs)
                self.fields['descriptor'].queryset = descriptorClass.objects.filter(calculatorSoftware='manual')
                self.fields['reaction'].widget = forms.HiddenInput()
                self.fields['reaction'].initial = reactionId
                self.fields['reaction'].queryset = Reaction.objects.filter(id=reactionId)
    
        return DescriptorValueForm

    return descriptorValueFormFactory
    
factory = descriptorValueFormFactoryFactory

NumRxnDescValFormFactory = factory(NumRxnDescriptorValue, NumRxnDescriptor)
OrdRxnDescValFormFactory = factory(OrdRxnDescriptorValue, OrdRxnDescriptor)
BoolRxnDescValFormFactory = factory(BoolRxnDescriptorValue, BoolRxnDescriptor)

def CatRxnDescValFormFactory(reactionId):
    """A specific factory for categorical descriptor value forms."""
    class CatRxnDescValForm(forms.ModelForm):

        """Special case because of the need to group permitted values."""
    
        value = forms.ChoiceField()
    
        class Meta:
            model = CatRxnDescriptorValue
            fields = ('descriptor', 'value', 'reaction')
    
        def __init__(self, *args, **kwargs):
            """Lock the reaction to the specified value."""
            super(CatRxnDescValForm, self).__init__(*args, **kwargs)
            descriptors = CatRxnDescriptor.objects.filter(calculatorSoftware='manual').prefetch_related('permittedValues')
            self.fields['descriptor'].queryset = descriptors
            self.fields['value'].choices = tuple(((descriptor.name, tuple((value.id, value.value) for value in descriptor.permittedValues.all())) for descriptor in descriptors))
            self.fields['reaction'].widget = forms.HiddenInput()
            self.fields['reaction'].initial = reactionId
            self.fields['reaction'].queryset = Reaction.objects.filter(id=reactionId)

    return CatRxnDescValForm
