'''A module containing forms for filtering compound objects'''
import django.forms as forms 
from django.core.exceptions import ValidationError
from django.forms.widgets import HiddenInput
from DRP.models import ChemicalClass



class CompoundFilterForm(forms.Form):
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

  def is_empty(self):
    return not any(self.cleaned_data.get(key) not in (None, '') for key in ('abbrev', 'CSID', 'INCHI', 'smiles')) or ('chemicalClasses' in self.cleaned_data and self.cleaned_data.get('chemicalClasses').count() != 0)

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

class CompoundFilterFormSet(forms.formsets.BaseFormSet):
  '''A formset for managing multiple filter forms, which OR together the results of each filter form to create a bigger queryset'''

  def __init__(self, user, labGroup, *args, **kwargs):
    self.form = CompoundFilterForm 
    self.user=user
    self.labGroup = labGroup
    self.extra=1
    self.can_delete = False
    self.can_order = False
    self.max_num=None
    self.validate_max=False
    self.absolute_max = 10000
    super(CompoundFilterFormSet, self).__init__(*args, **kwargs)

  def _construct_form(self, i, **kwargs):
    kwargs['user'] = self.user
    kwargs['labGroup'] = self.labGroup
    return super(CompoundFilterFormSet, self)._construct_form(i, **kwargs)

  @property
  def initData(self):
    initData = []
    for form in self:
      if form.is_valid() and not form.is_empty():
        initData.append(form.cleaned_data)
    return initData

  def fetch(self):
    qs = self.forms[0].fetch()
    for form in self:
      if not form.is_empty():
        qs |= form.fetch()
    return qs
