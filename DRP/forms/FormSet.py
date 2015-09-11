'''A formset class to make up for some of the deficiencies and standards
compliance issues raised by the django formset implementation'''

from widgets import SubmitButtonWidget
from django import forms
from django.core.exceptions import ValidationError
import copy

DELETE = "DELETEME"
TOTAL_FORMS = "TOTAL_FORMS"
MAX_FORMS = "MAX_FORMS"
PREFIX="PREFIX"
ADD_FORM = "ADD_FORM"
HARD_MAX_FORMS=10000

class FormSetManagerForm(forms.Form):

  def __init__(self, maxForms, formSetPrefix, initialCount, canAdd, canDelete, prefix='', data=None, *args, **kwargs):
    super(FormSetManagerForm, self).__init__(prefix=prefix+'-'+formSetPrefix, data=copy.copy(data), *args, **kwargs)
    self.fields[TOTAL_FORMS] = forms.IntegerField(min_value=0, widget=forms.widgets.HiddenInput, initial=1) #useful for javascript
    self.fields[MAX_FORMS] = forms.IntegerField(min_value=0, max_value=maxForms, initial=maxForms, widget=forms.widgets.HiddenInput) #useful for javascript
    self.fields[PREFIX] = forms.CharField(initial=formSetPrefix, widget=forms.widgets.HiddenInput, required=False) #useful for javascript
    if canAdd:
      self.canAdd = canAdd
      self.fields[ADD_FORM] = forms.BooleanField(label=None, widget=SubmitButtonWidget('Add another'), required=False)
    if canDelete:
      self.canDelete = canDelete
      self.fields[DELETE] = forms.BooleanField(label=None, widget=SubmitButtonWidget('Remove one'), required=False)

  def clean(self):
    cleaned_data = super(FormSetManagerForm, self).clean()
    if cleaned_data.get(TOTAL_FORMS) > cleaned_data.get(MAX_FORMS):
      raise ValidationError('No more may be added.', code='too_many_forms')
    if cleaned_data.get(ADD_FORM):
      cleaned_data[TOTAL_FORMS] +=1
    if cleaned_data.get(DELETE):
      cleaned_data[TOTAL_FORMS] -=1
    self.data[self.add_prefix(TOTAL_FORMS)] = str(self.cleaned_data.get(TOTAL_FORMS)) #YUK.
    if cleaned_data.get(TOTAL_FORMS) < 2 and (DELETE in self.fields):
      del self.fields[DELETE]
    self.cleaned_data = cleaned_data
    return cleaned_data

  @property
  def submittedForms(self):
    if self.cleaned_data.get(ADD_FORM):
      return self.cleaned_data.get(TOTAL_FORMS)-1
    else:
      return self.cleaned_data.get(TOTAL_FORMS)
      
  def as_ul(self):
    return self._html_output(
      normal_row = '<li%(html_class_attr)s>%(errors)s%(field)s%(help_text)s</li>',
      error_row = '<li>%s</li>',
      row_ender = "</li>",
      help_text_html = '<span calss="helptext">%s</span>',
      errors_on_separate_row=False
    )

class FormSet(object):

  def __init__(self, formClass, data=None, prefix='', initial=None, maxForms=HARD_MAX_FORMS, canDelete=False, canAdd=True):
    if initial is None:
      initialCount = 1
    elif len(initial) < 2:
      initialCount = 1
    else:
      initialCount = len(initial)
    self.managementForm = FormSetManagerForm(maxForms, prefix, initialCount, canAdd, canDelete, prefix='{}-manager'.format(prefix), data=data)
    if self.managementForm.is_bound and self.managementForm.is_valid():
        formCount = self.managementForm.submittedForms
    else:
      formCount = initialCount

    self.forms = []

    for i in range(0, formCount):
      self.forms.append(formClass(data=data, prefix='{}-{}'.format(prefix, i)))

    if self.managementForm.is_valid():
      if self.managementForm.cleaned_data.get(ADD_FORM):
        self.forms.append(formClass(prefix='{}-{}'.format(prefix, formCount)))

  def __iter__(self):
    for f in self.forms:
      yield f

  def __len__(self):
    return len(self.forms)

  def is_valid(self):
    return all(f.is_valid() for f in self.forms)

  @property
  def cleaned_data(self):
    return [f.cleaned_data for f in self.forms if f.is_valid()]

class ModelFormSet(FormSet):


  def __init__(self, modelClass, formClass=None, fields=None, *args, **kwargs):
    outerFields = fields
    if formClass is None:
      class ModelFormClass(forms.ModelForm):

        class Meta(object):
          model = modelClass
          if outerFields is not None:
            fields=outerFields
      formClass = ModelFormClass
    elif fields is not None:
      formClass._meta.fields=outerFields
    super(ModelFormSet, self).__init__(formClass, *args, **kwargs)

  def save(commmit=True):
    return [form.save(commit=commit) for form in self.forms if form.is_valid()]
