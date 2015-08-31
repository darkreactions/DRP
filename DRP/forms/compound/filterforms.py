'''A module containing forms for filtering compound objects'''
import django.forms as forms 
from django.core.exceptions import ValidationError
from django.forms.widgets import HiddenInput
from DRP.models import ChemicalClass, NumMolDescriptor, NumMolDescriptorValue, OrdMolDescriptor, OrdMolDescriptorValue
from DRP.models import CatMolDescriptor, CatMolDescriptorValue, BoolMolDescriptor, BoolMolDescriptorValue
from DRP.models import CatMolDescriptorPermitted
from DRP.forms import FilterForm, FilterFormSet, filterFormSetFactory

class CompoundFilterForm(FilterForm):
  '''A filter form to fetch Compound objects, a queryset of which is returned using the fetch() method.'''

  custom = forms.NullBooleanField(widget=forms.widgets.RadioSelect(choices=((None, 'Either'),(True, 'True'),(False, 'False'))), initial=None, required=False)

  def __init__(self, user, labGroup, *args, **kwargs):
    '''Sets up the form. Because most of the fields are based around models, they must be added dynamically.'''
    super(CompoundFilterForm, self).__init__(*args, **kwargs)
    self.empty_permitted = False #hard override to cope with a bad piece of programming in django.
    self.fields['abbrev'] = forms.ChoiceField(label='Abbreviation', choices=(('',''),) + tuple(((c['abbrev'],c['abbrev']) for c in labGroup.compound_set.all().values('abbrev').distinct())), required=False)
    self.fields['name'] = forms.ChoiceField(choices=((('',''),) + tuple((c['name'],c['name']) for c in labGroup.compound_set.all().values('name').distinct())), required=False)
    self.fields['chemicalClasses'] = forms.ModelMultipleChoiceField(label='Chemical Classes', queryset=ChemicalClass.objects.filter(compound__in=labGroup.compound_set.all()).distinct(), required=False)
    self.fields['CSID'] = forms.ChoiceField(choices=(('',''),) + tuple(((c['CSID'], c['CSID']) for c in labGroup.compound_set.all().values('CSID').distinct())), required=False)
    self.fields['INCHI'] = forms.CharField(required=False)
    self.fields['smiles'] = forms.CharField(required=False)
    self.fields['labGroup'] = forms.ModelChoiceField(queryset=user.labgroup_set.all(), initial=labGroup, widget=HiddenInput, error_messages={'invalid_choice':'You appear to have borrowed a search from a lab group to which you do not belong.'})
    self.fields['js_active'] = forms.NullBooleanField(widget=HiddenInput, required=False, initial=False)
    self.checkFields = ('abbrev', 'CSID', 'INCHI', 'smiles', 'chemicalClasses')

  def fetch(self):
    '''Fetches the labs according to data supplied. Expects the form to have been validated already.'''
    
    qs = self.cleaned_data['labGroup'].compound_set.all()
    if self.cleaned_data.get('js_active') not in ('', None, 'False'):
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
    if initial is None:
      init = {}
    else:
      init = initial #init points at the initial dictionary. What happens to init, happens to initial, and we need this for the pop methods.
    super(AdvancedCompoundFilterForm, self).__init__(initial=initial, *args, **kwargs) 
    self.numericFormSet = filterFormSetFactory(NumericFilterForm)(data=self.data, prefix = self.prefix+'_num', initial=init.pop('numeric', None))
    self.ordinalFormSet = filterFormSetFactory(OrdinalFilterForm)(data=self.data, prefix = self.prefix+'_ord', initial=init.pop('ordinal', None)),
    self.categoryFormSet = filterFormSetFactory(CategoryFilterForm)(data=self.data, prefix = self.prefix+'_cat', initial=init.pop('category', None)), 
    self.booleanFormSet = filterFormSetFactory(BooleanFilterForm)(data=self.data, prefix = self.prefix+'_bool', initial=init.pop('bool', None))]
    self.formSets = [self.numericFormSet, self.ordinalFormSet, self.categoryFormSet, self.booleanFormSet) 
                     
  def clean(self):
    cleaned_data = super(AdvancedCompoundFilterForm, self).clean()
    cleaned_data['numeric'] = self.numericFormSet.cleaned_data
    cleaned_data['ordinal'] = self.ordinalFormSet.cleaned_data
    cleaned_data['category'] = self.categoryFormSet.cleaned_data
    cleaned_data['booleanFormSet'] = self.booleanFormSet.cleaned_data

  def is_empty(self)
    empty = super(AdvancedCompoundFilterForm, self).is_empty()
    return all([empty] + [formSet.is_empty() for formSet in self.formSets]) 

  def is_valid(self):
    return super(AdvancedCompoundFilterForm, self).is_valid() and all(formSet.is_valid() for formSet in self.formSets) 

  def fetch(self):
    qs = super(AdvancedCompoundFilterForm, self).fetch()
    qs = qs.filter(nummoldescriptorvalue_set__in(self.numericFormSet.fetch()))
    qs = qs.filter(ordmoldescriptorvalue_set__in(self.ordinalFormSet.fetch()))
    qs = qs.filter(catmoldescriptorvalue_set__in(self.categoryFormSet.fetch())
    qs = qs.filter(boolmoldescriptorvalue_set__in(self.booleanFormSet.fetch())
    return qs

OPERATOR_CHOICES=(
  ('eq','='),
  ('gt','>'),
  ('ge', '&ge;'),
  ('lt', '<'),
  ('le', '&le;'),
  ('ne', '&ne;')
)

class QuantitativeFilterMixin:
  
  operator = forms.ChoiceField(choices=(OPERATOR_CHOICES))

  def applyFilters(self, qs):
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
      raise RuntimeError('Impossible Value provided to form, and passed validation')

class NumericFilterForm(FilterForm, QuantitativeFilterMixin):

  value = forms.DecimalField()

  def __init__(self, *args, **kwargs);
    super(NumericFilterForm, self).__init__(*args, **kwargs)
    self.fields['descriptor'] = forms.ModelChoiceField(queryset=NumMolDescriptor.objects.all())
    self.checkFields = ('value', 'descriptor')

  def fetch(self):
    qs = NumMolDescriptorValue.objects.filter(descriptor=self.cleaned_data.get('descriptor'))
    return self.applyFilters(qs)

def OrdinalFilterForm(FilterForm):

  def __init__(self, *args, **kwargs):
    super(OrdinalFilterForm, self).__init__(*args, **kwargs)
    self.fields['descriptor'] = forms.ModelChoiceField(queryset=OrdMolDescriptor.objects.all())
    self.fields['value'] = forms.ChoiceField(choices=((None, '') + tuple(md.heading, ((value, value) for value in range(md.minimum, md.maximum+1)) for md in OrdMolDescriptor.objects.all())))

  def fetch(self):
    qs = OrdMolDescriptorValue.objects.filter(self.cleaned_data.get('descriptor'))
    return self.applyFilters(qs)

def CategoryFilterForm(FilterForm)

  def __init__(self, *args, **kwargs):
    super(CategoryFilterForm, self).__init__(*args, **kwargs)
    self.fields['descriptor'] = forms.ModelChoiceField(queryset=CatMolDescriptor.objects.all())
    self.fields['value'] = forms.ChoiceField(choices=((None, '') + tuple(md.heading, ((value.pk, value.value) for value in CatMolDescriptorPermitted.objects.filter(descriptor=md)) for md in CatMolDescriptor.objects.all()))) #wow, that's hideous... It limits the options available to the available categorical molecular descriptor values, which are categorised according to the particular descriptor.

  def clean(self):
    super(CategoryFitlerForm, self).clean()
    if self.cleaned_data.get('value') not in self.cleaned_data.get('descriptor').permittedValues.all():
      raise ValidationError('The selected value is not appropriate to the selected descriptor', 'wrong_value')
    else:
      return self.cleaned_data

  def fetch(self):
    return CatMolDescriptorValue.objects.filter(descriptor=self.cleaned_data.get('descriptor'), value__pk==self.cleaned_data.get('value'))

def BooleanFilterForm(FilterForm):
  
  value = forms.NullBooleanField(widget=forms.widgets.RadioSelect(choices=((None, 'Either'),(True, 'True'),(False, 'False'))), initial=None, required=False)
  
  def __init__(self, *args, **kwargs):
    super(BooleanFilterForm, self).__init__(*args, **kwargs)
    self.fields['descriptor'] = forms.ModelChoiceField(queryset=BoolMolDescriptor.objects.all())

  def fetch(self):
    return BoolMolDescriptorValue.objects.filter(descriptor=self.cleaned_data.get('descriptor'), value=self.cleaned_data.get('value'))

class CompoundFilterFormSet(FilterFormSet):
  '''A formset for managing multiple filter forms, which OR together the results of each filter form to create a bigger queryset'''

  form = CompoundFilterForm 

  def __init__(self, user, labGroup, *args, **kwargs):
    self.user=user
    self.labGroup = labGroup
    super(CompoundFilterFormSet, self).__init__(*args, **kwargs)

  def _construct_form(self, i, **kwargs):
    kwargs['user'] = self.user
    kwargs['labGroup'] = self.labGroup
    return super(CompoundFilterFormSet, self)._construct_form(i, **kwargs)

class AdvancedCompoundFilterFormSet(CompoundFilterFormSet):

  form = AdvancedCompoundFilterForm
