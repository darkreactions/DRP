'''Forms for the django admin for custom descriptors'''

from django import forms
from django.core.exceptions import ValidationError
from DRP.models import CatRxnDescriptor, OrdRxnDescriptor, NumRxnDescriptor, BoolRxnDescriptor
from DRP.models import CategoricalDescriptorPermittedValue, CategoricalDescriptor


class DescriptorAdmin(forms.ModelForm):
    '''A mixin for behaviours common to all descriptor admin forms'''

    def clean(self, *args, **kwargs):
        '''This clean method is purely desingned to stop the overwriting of plugin descriptors'''
        cleaned_data = super(DescriptorAdmin, self).clean(*args, **kwargs)
        if self.instance.pk and self.instance.calculatorSoftware != 'manual':
            raise ValidationError('This descriptor is not a manual descriptor, and thus cannot be edited using the django admin', 'not_manual')
        return cleaned_data

    def save(self, commit=True, *args, **kwargs):
        descriptor = super(DescriptorAdmin, self).save(commit=False, *args, **kwargs)
        descriptor.calculatorSoftware = 'manual'
        descriptor.calculatorSoftwareVersion = '0'
        if commit:
            descriptor.save()
        return descriptor


class CatRxnDescriptorForm(DescriptorAdmin):
    '''An admin form for custom Categorical Reaction Descriptors'''

    class Meta:
        fields = ('heading', 'name')
        model = CatRxnDescriptor


class CatDescPermittedValueForm(forms.ModelForm):
    '''A mechanism to create permitted values for custom Categorical Reaction descriptors'''

    class Meta:
        model = CategoricalDescriptorPermittedValue
        fields = ('descriptor', 'value')

    def clean(self, *args, **kwargs):
        cleaned_data = super(CatDescPermittedValueForm, self).clean(*args, **kwargs)
        if not CategoricalDescriptor.objects.filter(calculatorSoftware='manual', descriptor=self.instance.descriptor).exists():
            raise ValidationError('You may only edit descriptor values for your own custom descriptors')
        else:
            return cleaned_data

    def __init__(self, *args, **kwargs):
        super(CatDescPermittedValueForm, self).__init__(*args, **kwargs)
        self.fields['descriptor'].queryset = CategoricalDescriptor.objects.filter(calculatorSoftware='manual')


class OrdRxnDescriptorForm(DescriptorAdmin):
    '''An admin form for creating custom Ordinal reaction descriptors'''

    class Meta:
        fields = ('heading', 'name', 'minimum', 'maximum')
        model = OrdRxnDescriptor


class NumRxnDescriptorForm(DescriptorAdmin):
    '''An admin form for creating custom numeric reaction descriptors'''

    class Meta:
        fields = ('heading', 'name', 'minimum', 'maximum')
        model = NumRxnDescriptor


class BoolRxnDescriptorForm(DescriptorAdmin):
    '''An admin form for creating custom boolean reaction descriptors'''

    class Meta:
        fields = ('heading', 'name')
        model = BoolRxnDescriptor
