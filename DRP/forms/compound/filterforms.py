'''A module containing forms for filtering compound objects'''
import django.forms as forms 
from django.conf import settings
from django.core.exceptions import ValidationError
from django.forms.widgets import HiddenInput
from DRP.models import Compound, ChemicalClass, NumMolDescriptor, NumMolDescriptorValue, OrdMolDescriptor, OrdMolDescriptorValue
from DRP.models import CatMolDescriptor, CatMolDescriptorValue, BoolMolDescriptor, BoolMolDescriptorValue
from DRP.models import CatMolDescriptorPermitted
from DRP.forms import FilterForm, FilterFormSet, filterFormSetFactory
from django.utils.safestring import mark_safe
from operator import and_

class CompoundFilterForm(FilterForm):
  '''A filter form to fetch Compound objects, a queryset of which is returned using the fetch() method.'''

  model = Compound
  custom = forms.NullBooleanField(widget=forms.widgets.RadioSelect(choices=((None, 'Either'),(True, 'True'),(False, 'False'))), initial=None, required=False)

  def __init__(self, user, labGroup, *args, **kwargs):
    '''Sets up the form. Because most of the fields are based around models, they must be added dynamically.'''
    super(CompoundFilterForm, self).__init__(*args, **kwargs)
    self.empty_permitted = False #hard override to cope with a bad piece of programming in django.
    self.fields['abbrev'] = forms.ChoiceField(label='Abbreviation', choices=(('',settings.EMPTY_LABEL),) + tuple(((c['abbrev'],c['abbrev']) for c in labGroup.compound_set.all().values('abbrev').distinct())), required=False)
    self.fields['name'] = forms.ChoiceField(choices=((('', settings.EMPTY_LABEL),) + tuple((c['name'],c['name']) for c in labGroup.compound_set.all().values('name').distinct())), required=False)
    self.fields['chemicalClasses'] = forms.ModelMultipleChoiceField(label='Chemical Classes', queryset=ChemicalClass.objects.filter(compound__in=labGroup.compound_set.all()).distinct(), required=False)
    self.fields['CSID'] = forms.ChoiceField(choices=(('',settings.EMPTY_LABEL),) + tuple(((c['CSID'], c['CSID']) for c in labGroup.compound_set.all().values('CSID').distinct())), required=False)
    self.fields['INCHI'] = forms.CharField(required=False)
    self.fields['smiles'] = forms.CharField(required=False)
    self.fields['labGroup'] = forms.ModelChoiceField(queryset=user.labgroup_set.all(), initial=labGroup, widget=HiddenInput, error_messages={'invalid_choice':'You appear to have borrowed a search from a lab group to which you do not belong.'}, empty_label=settings.EMPTY_LABEL)
    self.fields['js_active'] = forms.NullBooleanField(widget=HiddenInput, required=False, initial=False)
    self.checkFields = ('abbrev', 'CSID', 'INCHI', 'smiles')

  def is_empty(self):
    '''Checks that the form is empty and performs specific checks for chemicalClasses'''
    base_empty = super(CompoundFilterForm, self).is_empty() #performs the normal check on the easy fields
    return base_empty and (self.cleaned_data.get('chemicalClasses') not in self.fields['chemicalClasses'].empty_values or self.cleaned_data.get('chemicalClasses').count() == 0)

  def fetch(self):
    '''Fetches the labs according to data supplied. Exp the form to have been validated already.'''
    
    qs = self.cleaned_data['labGroup'].compound_set.all()
    if self.cleaned_data.get('js_active') not in ('', None, False):
      raise RuntimeError(self.cleaned_data.get('js_active'))
    else:
      if self.cleaned_data.get('abbrev') not in (None, ''):
        qs = qs.filter(abbrev__contains=self.cleaned_data['abbrev'])
      if self.cleaned_data.get('name') not in (None, ''):
        qs = qs.filter(name__contains=self.cleaned_data['name'])
    if self.cleaned_data['chemicalClasses'].count() != 0:
      for cc in self.cleaned_data['chemicalClasses']:
        qs = qs.filter(chemicalClasses=cc)
    if self.cleaned_data.get('CSID') not in (None, ''):
      qs = qs.filter(CSID=self.cleaned_data['CSID'])
    if self.cleaned_data.get('INCHI') not in (None, ''):
      qs = qs.filter(INCHI=self.cleaned_data['INCHI'])
    if self.cleaned_data.get('smiles') not in (None, ''):
      qs = qs.filter(smiles=self.cleaned_data['smiles'])
    if self.cleaned_data.get('custom') not in (None, ''):
      qs = qs.filter(custom=True if self.cleaned_data.get('custom') is 'True' else False)
    return qs

class AdvancedCompoundFilterForm(CompoundFilterForm):
  '''A form for making more complex queries about compounds, specifically using their descriptor values'''

  def __init__(self, initial=None, *args, **kwargs):
    '''Sets up FormSets for managing the descriptor filters as a part of this form'''
    if initial is None:
      init = {}
    else:
      init = initial #init points at the initial dictionary. What happens to init, happens to initial, and we need this for the pop methods.
    super(AdvancedCompoundFilterForm, self).__init__(initial=initial, *args, **kwargs) 
    data = None if self.data=={} else self.data
    self.numericFormSet = filterFormSetFactory(NumericFilterForm, NumMolDescriptorValue)(data=data, prefix = '{}_num'.format(self.prefix), initial=init.pop('numeric', None), operator=and_)
    self.ordinalFormSet = filterFormSetFactory(OrdinalFilterForm, OrdMolDescriptorValue)(data=data, prefix = '{}_ord'.format(self.prefix), initial=init.pop('ordinal', None), operator=and_)
    self.categoryFormSet = filterFormSetFactory(CategoryFilterForm, CatMolDescriptorValue)(data=data, prefix = '{}_cat'.format(self.prefix), initial=init.pop('category', None), operator=and_)
    self.booleanFormSet = filterFormSetFactory(BooleanFilterForm, BoolMolDescriptorValue)(data=data, prefix = '{}_bool'.format(self.prefix), initial=init.pop('bool', None), operator=and_)
    self.formSets = [self.numericFormSet, self.ordinalFormSet, self.categoryFormSet, self.booleanFormSet] 
                     
  def clean(self):
    '''Returns cleaned data, with additional data appended for this forms attached formsets'''
    cleaned_data = super(AdvancedCompoundFilterForm, self).clean()
    cleaned_data['numeric'] = self.numericFormSet.cleaned_data
    cleaned_data['ordinal'] = self.ordinalFormSet.cleaned_data
    cleaned_data['category'] = self.categoryFormSet.cleaned_data
    cleaned_data['booleanFormSet'] = self.booleanFormSet.cleaned_data
    return cleaned_data

  def is_empty(self):
    '''Checks that the form is empty and performs special checks for the attached formsets'''
    empty = super(AdvancedCompoundFilterForm, self).is_empty()
    return all([empty] + [formSet.is_empty() for formSet in self.formSets]) 

  def is_valid(self):
    '''validates the formsets as well as this form'''
    #raise RuntimeError([(formSet.is_valid(), formSet._errors) for formSet in self.formSets])
    return super(AdvancedCompoundFilterForm, self).is_valid() and all(formSet.is_valid() for formSet in self.formSets) 

  def fetch(self):
    '''returns compounds as per the filters'''
    qs = super(AdvancedCompoundFilterForm, self).fetch()
    qs = qs.filter(nummoldescriptorvalue__in=self.numericFormSet.fetch())
    qs = qs.filter(ordmoldescriptorvalue__in=self.ordinalFormSet.fetch())
    qs = qs.filter(catmoldescriptorvalue__in=self.categoryFormSet.fetch())
    qs = qs.filter(boolmoldescriptorvalue__in=self.booleanFormSet.fetch())
    return qs

OPERATOR_CHOICES=(
  ('eq','='),
  ('gt','>'),
  ('ge', mark_safe('&ge;')),
  ('lt', '<'),
  ('le', mark_safe('&le;')),
  ('ne', mark_safe('&ne;'))
)

class QuantitativeFilterMixin(forms.Form):
  '''A mixin which contains information which is used by both NumericFilterForms and OrdinalFilterForms'''
  
  operator = forms.ChoiceField(choices=(OPERATOR_CHOICES))

  def applyFilters(self, qs):
    '''Works applies the correct filtration operator to the queryset'''
    op = self.cleaned_data.get('operator')
    value = self.cleaned_data.get('value')
    if op == 'eq':
      return qs.filter(value=value)
    elif op == 'gt':
      return qs.filter(value__gt=value)
    elif op == 'ge':
      return qs.filter(value__gte=value)
    elif op == 'lt':
      return qs.filter(value__lt=value)
    elif op == 'le':
      return qs.filter(value__lte=value)
    elif op == 'ne':
      return qs.filter(value__ne=value)
    else:
      raise RuntimeError('Impossible Value provided to form, and passed validation:{}'.format(op))

class NumericFilterForm(QuantitativeFilterMixin, FilterForm):
  '''A form to obtain numeric descriptor values. Used by the advanced compound filtering form'''

  def __init__(self, *args, **kwargs):
    '''Sets up the forms fields, almost all which require some level of dynamism'''
    super(NumericFilterForm, self).__init__(*args, **kwargs)
    self.fields['descriptor'] = forms.ModelChoiceField(queryset=NumMolDescriptor.objects.all(), required=False, empty_label=settings.EMPTY_LABEL)
    self.fields['value'] = forms.DecimalField(required=False)
    self.fields['operator'] = self.fields.pop('operator') #because there isn't a more sensible way of doing the re-ordering for non model forms.
    self.checkFields = ('value', 'descriptor')

  def clean(self):
    '''checks that the descriptor choice and value have both been supplied, or neither have'''
    super(NumericFilterForm, self).clean()
    if (self.cleaned_data.get('descriptor') is None) ^ (self.cleaned_data.get('value') is None):
      raise ValidationError('Both a descriptor and a value must be provided. Empty the fields completely to ignore this input.')
    else:
      return self.cleaned_data

  def fetch(self):
    '''fetch the NumMolDescriptor objects which match the supplied form values'''
    qs = NumMolDescriptorValue.objects.filter(descriptor=self.cleaned_data.get('descriptor'))
    return self.applyFilters(qs)

  def is_empty(self):
    '''returns true if all form fields were empty at submission'''
    empty = super(NumericFilterForm, self).is_empty()
    return empty and self.cleaned_data.get('descriptor') is None

class OrdinalFilterForm(QuantitativeFilterMixin, FilterForm):
  '''A form to obtain Ordinal Descriptor values, used by the advanced compound filtering form'''

  def __init__(self, *args, **kwargs):
    '''Sets the fields up in the right order'''
    super(OrdinalFilterForm, self).__init__(*args, **kwargs)
    self.fields['descriptor'] = forms.ModelChoiceField(queryset=OrdMolDescriptor.objects.all(), required=False, empty_label=settings.EMPTY_LABEL)
    self.fields['value'] = forms.ChoiceField(choices=((('', settings.EMPTY_LABEL),) + tuple((md.name, tuple((value, value) for value in range(md.minimum, md.maximum+1))) for md in OrdMolDescriptor.objects.all())), required=False)
    self.fields['operator'] = self.fields.pop('operator') #because there isn't a more sensible way of doing the re-ordering for non model forms.
    self.checkFields = ('value', 'descriptor')

  def clean(self):
    super(OrdinalFilterForm, self).clean()
    '''checks that the descriptor choice and value have both been supplied, or neither have'''
    if (self.cleaned_data.get('descriptor') is None) ^ (self.cleaned_data.get('value') == ''):
      raise ValidationError('Both a descriptor and a value must be provided. Empty the fields completely to ignore this input.')
    else:
      return self.cleaned_data

  def fetch(self):
    '''fetch the NumMolDescriptor objects which match the supplied form values'''
    qs = OrdMolDescriptorValue.objects.filter(self.cleaned_data.get('descriptor'))
    return self.applyFilters(qs)

  def is_empty(self):
    '''returns true if all form fields were empty at submission'''
    empty = super(OrdinalFilterForm, self).is_empty()
    return empty and self.cleaned_data.get('descriptor') is None

class CategoryFilterForm(FilterForm):
  '''A filter form for obtaining Categorical descriptor values. Used by the advanced compound filtering form'''

  def __init__(self, *args, **kwargs):
    '''Sets teh forms up in the right order'''
    super(CategoryFilterForm, self).__init__(*args, **kwargs)
    self.fields['descriptor'] = forms.ModelChoiceField(queryset=CatMolDescriptor.objects.all(), required=False, empty_label=settings.EMPTY_LABEL)
    self.fields['value'] = forms.ChoiceField(choices=((('', settings.EMPTY_LABEL),) + tuple((md.name, tuple((value.pk, value.value) for value in CatMolDescriptorPermitted.objects.filter(descriptor=md))) for md in CatMolDescriptor.objects.all())), required=False) #wow, that's hideous... It limits the options available to the available categorical molecular descriptor values, which are categorised according to the particular descriptor.
    self.checkFields = ('value', 'descriptor')

  def clean(self):
    '''Checks that both or neither of value and descriptor have been supplied, and checks that a descriptor value choice appropriate to the
    descriptor has been chosen'''
    super(CategoryFilterForm, self).clean()
    if (self.cleaned_data.get('descriptor') is None) ^ (self.cleaned_data.get('value') is ''):
      raise ValidationError('Both a descriptor and a value must be provided. Empty the fields completely to ignore this input.')
    if self.cleaned_data.get('descriptor') is not None and CatMolDescriptorValue.objects.get(pk=self.cleaned_data.get('value')) not in self.cleaned_data.get('descriptor').permittedValues.all():
      raise ValidationError('The selected value is not appropriate to the selected descriptor', 'wrong_value')
    else:
      return self.cleaned_data

  def fetch(self):
    '''returns the categorical descriptor value objects'''
    return CatMolDescriptorValue.objects.filter(descriptor=self.cleaned_data.get('descriptor'), value__pk=self.cleaned_data.get('value'))

  def is_empty(self):
    '''returns true if all form fields were empty at submission'''
    empty = super(CategoryFilterForm, self).is_empty()
    return empty and self.cleaned_data.get('descriptor') is None and self.cleaned_data.get('value') is None

class BooleanFilterForm(FilterForm):
  '''A form for filtering boolean descriptor values.'''
  
  def __init__(self, *args, **kwargs):
    '''Sets teh fields up for this form'''
    super(BooleanFilterForm, self).__init__(*args, **kwargs)
    self.fields['descriptor'] = forms.ModelChoiceField(queryset=BoolMolDescriptor.objects.all(), required=False, empty_label=settings.EMPTY_LABEL)
    self.fields['value'] = forms.NullBooleanField(widget=forms.widgets.RadioSelect(choices=((None, 'Either'),(True, 'True'),(False, 'False'))), initial=None, required=False)
    self.checkFields = ('value', 'descriptor')

  def fetch(self):
    '''returns the appropriate queryset'''
    return BoolMolDescriptorValue.objects.filter(descriptor=self.cleaned_data.get('descriptor'), value=self.cleaned_data.get('value'))

  def is_empty(self):
    '''returns true if all form fields were empty at submission'''
    empty = super(BooleanFilterForm, self).is_empty()
    return empty and self.cleaned_data.get('value') is None and self.cleaned_data.get('descriptor') is None

class CompoundFilterFormSet(FilterFormSet):
  '''A formset for managing multiple filter forms, which OR together the results of each filter form to create a bigger queryset'''

  form = CompoundFilterForm 

  def __init__(self, user, labGroup, *args, **kwargs):
    '''Initialises the formset with the user and lab group variables needed to construct the forms'''
    self.user=user
    self.labGroup = labGroup
    super(CompoundFilterFormSet, self).__init__(*args, **kwargs)

  def _construct_form(self, i, **kwargs):
    '''constructs the forms using an overriden private method'''
    kwargs['user'] = self.user
    kwargs['labGroup'] = self.labGroup
    return super(CompoundFilterFormSet, self)._construct_form(i, **kwargs)

class AdvancedCompoundFilterFormSet(CompoundFilterFormSet):
  '''A formset for advanced compound filtering'''

  form = AdvancedCompoundFilterForm
