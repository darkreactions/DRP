'''A module containing forms for filtering compound objects'''
import django.forms as forms 
from django.conf import settings
from django.core.exceptions import ValidationError
from django.forms.widgets import HiddenInput
from DRP.models import Compound, ChemicalClass, NumMolDescriptor, NumMolDescriptorValue, OrdMolDescriptor, OrdMolDescriptorValue
from DRP.models import CatMolDescriptor, CatMolDescriptorValue, BoolMolDescriptor, BoolMolDescriptorValue
from DRP.models import CategoricalDescriptorPermittedValue
from DRP.forms import FilterForm, FilterFormSet, filterFormSetFactory
from django.utils.safestring import mark_safe
from operator import and_

OPERATOR_CHOICES=(
  ('eq','='),
  ('gt','>'),
  ('ge', mark_safe('&ge;')),
  ('lt', '<'),
  ('le', mark_safe('&le;')),
  ('ne', mark_safe('&ne;'))
)

class PerformedReactionFilterForm(FilterForm):
  '''A filter form to fetch Compound objects, a queryset of which is returned using the fetch() method.'''

  model = PerformedReaction
  custom = forms.NullBooleanField(widget=forms.widgets.RadioSelect(choices=((None, 'Either'),(True, 'True'),(False, 'False'))), initial=None, required=False)

  def __init__(self, user, labGroup, *args, **kwargs):
    '''Sets up the form. Because most of the fields are based around models, they must be added dynamically.'''
    super(PerformedReactionsFilterForm, self).__init__(*args, **kwargs)
    self.empty_permitted = False #hard override to cope with a bad piece of programming in django.
    self.fields['reference']=forms.CharField(required=False)
    self.fields['legacyRecommendedFlag'] = forms.NullBooleanField(widget=forms.widgets.RadioSelect(choices((None, 'Either'), (True, 'True'), (False, 'False'))), initial=None, required=false) 
    self.fields['valid'] = forms.NullBooleanField(widget=forms.widgets.RadioSelect(choices((None, 'Either'), (True, 'True'), (False, 'False'))), initial=None, required=False) 
    self.fields['public'] = forms.NullBooleanField(widget=forms.widgets.RadioSelect(choices((None, 'Either'), (True, 'True'), (False, 'False'))), initial=None, required=False) 
    self.fields['duplicateOf'] = forms.modelChoiceField(queryset=PerformedReactions.objects.filter(labgroup__in=user.labgroup_set), required=False)
    self.fields['labGroup'] = forms.ModelChoiceField(queryset=user.labgroup_set.all(), initial=labGroup, widget=HiddenInput, error_messages={'invalid_choice':'You appear to have borrowed a search from a lab group to which you do not belong.'}, empty_label=settings.EMPTY_LABEL)
    self.fields['js_active'] = forms.NullBooleanField(widget=HiddenInput, required=False, initial=False)
    self.checkFields = ('reference', 'recommendation', 'legacyRecommendedFlag', 'valid', 'public', 'duplicateOf', 'inTrainingSetFor', 'inTestSetFor','labGroup')

  def is_empty(self):
    '''Checks that the form is empty and performs specific checks for chemicalClasses'''
    base_empty = super(PerformedReactionFilterForm, self).is_empty() #performs the normal check on the easy fields
    return base_empty 

  def fetch(self):
    '''Fetches the compounds according to data supplied. Exp the form to have been validated already.'''
    
    qs = self.cleaned_data['labGroup'].PerformedReaction.objects.filter(labgroup__in=user.labgroup_set)
    if self.cleaned_data.get('js_active') not in ('', None, False):
      raise RuntimeError(self.cleaned_data.get('js_active'))
    else:
      if self.cleaned_data.get('reference') not in (None, ''):
        qs = qs.filter(reference=self.cleaned_data['reference'])
      if self.cleaned_data.get('recommendation') not in (None, ''):
        qs = qs.filter(recommendation=self.cleaned_data['recommendation'])
      if self.cleaned_data.get('legacyRecommendedFlag') not in (None, ''):
        qs = qs.filter(legacyRecommendedFlag=self.cleaned_data['legacyRecommendedFlag'])
      if self.cleaned_data.get('valid') not in (None, ''):
        qs = qs.filter(valid=self.cleaned_data['valid'])
      if self.cleaned_data.get('public') not in (None, ''):
        qs = qs.filter(public=self.cleaned_data['public'])
      if self.cleaned_data.get('duplicateOf') not in (None, ''):
        qs = qs.filter(duplicateOf=self.cleaned_data['duplicateOf'])
      if self.cleaned_data.get('custom') not in (None, ''):
        qs = qs.filter(custom=True if self.cleaned_data.get('custom') is 'True' else False)
      return qs
  
class AdvancedPerformedReactionFilterForm(PerformedReactionFilterForm):
  '''A form for making more complex queries about compounds, specifically using their descriptor values'''

  def __init__(self, initial=None, *args, **kwargs):
    '''Sets up FormSets for managing the descriptor filters as a part of this form'''
    if initial is None:
      init = {}
    else:
      init = initial #init points at the initial dictionary. What happens to init, happens to initial, and we need this for the pop methods.
    super(AdvancedPerformedReactionFilterForm, self).__init__(initial=initial, *args, **kwargs) 
    data = None if self.data=={} else self.data
    self.numericFormSet = filterFormSetFactory(NumericFilterForm, NumMolDescriptorValue)(data=data, prefix = '{}_num'.format(self.prefix), initial=init.pop('numeric', None), operator=and_)
    self.ordinalFormSet = filterFormSetFactory(OrdinalFilterForm, OrdMolDescriptorValue)(data=data, prefix = '{}_ord'.format(self.prefix), initial=init.pop('ordinal', None), operator=and_)
    self.categoryFormSet = filterFormSetFactory(CategoryFilterForm, CatMolDescriptorValue)(data=data, prefix = '{}_cat'.format(self.prefix), initial=init.pop('category', None), operator=and_)
    self.booleanFormSet = filterFormSetFactory(BooleanFilterForm, BoolMolDescriptorValue)(data=data, prefix = '{}_bool'.format(self.prefix), initial=init.pop('bool', None), operator=and_)
    self.compoundQuantityFormSet = filterFormSetFactory(CompoundQuantityFilterForm, CompoundQuantity)(data=data, prefix = '{}_compoundquantity'.format(self.prefix), initial=init.pop('compoundquantity', None), operator=and_)
    self.performedDateFormSet = filterFormSetFactory(PerformedDateFilterForm, PerformedReaction)(date=date,prefix='{}_performeddate'.format(self.prefix), initial=init.pop('performeddate'), None), operator=and_)
    self.insertedDateFormSet = filterFormSetFactory(InsertedDateFilterForm, PerformedReaction)(date=date, prefix='{}_inserteddate'.format(self.prefix), initial=init.pop('inserteddate'), None), operator=and_)
    self.formSets = [self.numericFormSet, self.ordinalFormSet, self.categoryFormSet, self.booleanFormSet, self.compoundquantityFormSet] 
    
  def clean(self):
    '''Returns cleaned data, with additional data appended for this forms attached formsets'''
    cleaned_data = super(AdvancedPerformedReactionFilterForm, self).clean()
    cleaned_data['numeric'] = self.numericFormSet.cleaned_data
    cleaned_data['ordinal'] = self.ordinalFormSet.cleaned_data
    cleaned_data['category'] = self.categoryFormSet.cleaned_data
    cleaned_data['booleanFormSet'] = self.booleanFormSet.cleaned_data
    cleaned_data['compoundquantity'] = self.compoundquantityFormSet.cleaned_data
    return cleaned_data

  def is_empty(self):
    '''Checks that the form is empty and performs special checks for the attached formsets'''
    empty = super(AdvancedPerformedReactionFilterForm, self).is_empty()
    return all([empty] + [formSet.is_empty() for formSet in self.formSets]) 

  def is_valid(self):
    '''validates the formsets as well as this form'''
    #raise RuntimeError([(formSet.is_valid(), formSet._errors) for formSet in self.formSets])
    return super(AdvancedPerformedReactionFilterForm, self).is_valid() and all(formSet.is_valid() for formSet in self.formSets) 

  def fetch(self):
    '''returns compounds as per the filters'''
    qs = super(AdvancedPerformedReactionFilterForm, self).fetch()
    qs = qs.filter(numrxndescriptorvalue__in=self.numericFormSet.fetch())
    qs = qs.filter(ordrxndescriptorvalue__in=self.ordinalFormSet.fetch())
    qs = qs.filter(catrxndescriptorvalue__in=self.categoryFormSet.fetch())
    qs = qs.filter(boolrxndescriptorvalue__in=self.booleanFormSet.fetch())
    qs = qs.filter(compoundquantity__in=self.compoundQuantityFormSet.fetch())
    qs = qs & performedDateFormSet.fetch()
    qs = qs & insertedDateFormSet.fetch() 
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
    # if False ^ False (meaning that both are supplied), skip to the else (no error)
    # if True ^ True, then False (according to xor logic), skip to else (no error)
    # Otherwise, if False ^ True | if True ^ False, then if True --> raise error 
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
    qs = OrdMolDescriptorValue.objects.filter(descriptor=self.cleaned_data.get('descriptor'))
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
    self.fields['value'] = forms.ChoiceField(choices=((('', settings.EMPTY_LABEL),) + tuple((md.name, tuple((value.pk, value.value) for value in CategoricalDescriptorPermittedValue.objects.filter(descriptor=md))) for md in CatMolDescriptor.objects.all())), required=False) #wow, that's hideous... It limits the options available to the available categorical molecular descriptor values, which are categorised according to the particular descriptor.
    self.checkFields = ('value', 'descriptor')

  def clean(self):
    '''Checks that both or neither of value and descriptor have been supplied, and checks that a descriptor value choice appropriate to the
    descriptor has been chosen'''
    super(CategoryFilterForm, self).clean()
    if (self.cleaned_data.get('descriptor') is None) ^ (self.cleaned_data.get('value') is ''):
      raise ValidationError('Both a descriptor and a value must be provided. Empty the fields completely to ignore this input.')
    if self.cleaned_data.get('descriptor') is not None and self.cleaned_data.get('descriptor').permittedValues.filter(pk=self.cleaned_data.get("value")).count()<1:
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


class CompoundQuantityFilterForm(FilterForm):
  '''A formset for filtering compound quantities, roles, and amounts.'''
  def __init__(self, user, *args, **kwargs):
    super(CompoundQuantityFilterForm, self).__init__(*args, **kwargs)
    self.fields['compound'] = forms.ModelChoiceField(queryset=CompoundQuantity.objects.filter(compound__labgroup__in = user.labgroup_set), required=False, empty_label=settings.EMPTY_LABEL)
    self.fields['role'] = forms.ModelChoiceField(queryset=CompoundQuantity.objects.all(), required=False, empty_label=settings.EMPTY_LABEL)
    self.fields['amount'] = forms.FloatField(required=False, empty_label=settings.EMPTY_LABEL)

  def fetch(self):
    '''returns the appropriate queryset'''
    return CompoundQuantity.objects.filter(compound=self.cleaned_data.get('compound'), role=self.cleaned_data.get('role'), amount=self.cleaned_data.get('amount')) 

  def is_empty(self):
    '''returns true if all form fields were empty at submission'''
    empty = super(CompoundQuantityFilterForm, self).is_empty()
    return empty and self.cleaned_data.get('compound') is None and self.cleaned_data.get('role') is None and slef.cleaned_data.get('amount') is None 


DATE_OPERATOR_CHOICES = (
  ('eq','on'),
  ('gt','after'),
  ('ge', 'on or after'),
  ('lt', 'before'),
  ('le', 'on or before'),
  ('ne', 'everything other than'))

class PerformedDateFilterForm(QuantitativeMixin, FilterForm): 
  '''A formset for filtering reactions by date performed'''
  operator = forms.ChoiceField(choices=(DATE_OPERATOR_CHOICES))
  date = forms.DateField(required=False)
  def __init__(self, user, *args, **kwargs):
    super(PerformedDateFilterForm).__init__(*args, **kwargs)
    self.user=user

  def clean(self):
    '''checks that the descriptor choice and value have both been supplied, or neither have'''
    super(PerformedDateFilterForm, self).clean()
    # if False ^ False (meaning that both are supplied), skip to the else (no error)
    # if True ^ True, then False (according to xor logic), skip to else (no error)
    # Otherwise, if False ^ True | if True ^ False, then if True --> raise error 
    if (self.cleaned_data.get('date') is None) ^ (self.cleaned_data.get('operator') is not None):
      raise ValidationError('Both a descriptor and a value must be provided. Empty the fields completely to ignore this input.')
    else:
      return self.cleaned_data
   
  def fetch()
   qs = PerformedReaction.objects.filter(performeddate=self.cleaned_data.get('date'))
    return self.applyFilters(qs)

  def is_empty()
    '''returns true if all form fields were empty at submission'''
    empty = super(PerformedDateFilterForm, self).is_empty()
    return empty and self.cleaned_data.get('date') is None

class InsertedDateFilterForm(PerformedDateFilterForm):
  '''A formset for filtering reactions by date inserted'''
 
  def fetch()
   qs = PerformedReaction.objects.filter(inserteddate=self.cleaned_data.get('date'))
    return self.applyFilters(qs)





class PerformedReactionFilterFormSet(FilterFormSet):
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

class AdvancedPerformedReactionFilterFormSet(PerformedReactionFilterFormSet):
  '''A formset for advanced compound filtering'''

  form = AdvancedCompoundFilterForm
