'''A formset class to make up for some of the deficiencies and standards
compliance issues raised by the django formset implementation'''

from widgets import SubmitButtonWidget
from django import forms
from django.core.exceptions import ValidationError

DELETE = "DELETEME"
TOTAL_FORMS = "TOTAL_FORMS"
MAX_FORMS = "MAX_FORMS"
PREFIX="PREFIX"
ADD_FORM = "ADD_FORM"
HARD_MAX_FORMS=10000

class FormSetManagerForm(forms.Form)

  def __init__(self, maxForms, formSetPrefix, initialCount, canAdd, canDelete, prefix='', *args, **kwargs)
    super(FormSetManagerForm, self).__init__(prefix=prefix+'-'+formSetPrefix, *args, **kwargs)
    self.fields[TOTAL_FORMS] = forms.IntegerField(min_value=0, widget=forms.widgets.HiddenInput, initial=) #useful for javascript
    self.fields[MAX_FORMS] = forms.IntegerField(min_value=0, max_value=maxForms, initial=maxForms) #useful for javascript
    self.fields[PREFIX] = forms.CharField(initial=formset_prefix) #useful for javascript
    if canAdd:
      self.canAdd = canAdd
      self.fields[ADD_FORM] = forms.BooleanField(widget=SubmitButtonWidget(value='Add Another'))

    def clean(self):
      cleaned_data = super(FormSetManagerForm, self).clean()
      if TOTAL_FORMS in self.cleaned_data:
          cleaned_data[TOTAL_FORMS] += 1
        if cleaned_data.get(TOTAL_FORMS) > cleaned_data.get(MAX_FORMS):
          raise ValidationError('No more may be added.', code='too_many_forms')

class FormSet(object):

  def __init__(self, FormClass, data=None, prefix='', initial=None, maxForms=HARD_MAX_FORMS, canDelete=False, canAdd=True):
    if initial is None:
      initialCount = 1
    elif len(initial) < 2:
      initialCount = 1
    else:
      initialCount = len(initial)
    self.managementForm = FormSetManagerForm(maxForms, prefix, initialCount, canAdd, canDelete, prefix='{}-manager'.format(prefix), data=data)
    if self.managementForm.is_bound() and self.managementForm.is_valid():
        formCount = self.managementForm.cleaned_data.get(TOTAL_FORMS) - (1 if self.managementForm.get(ADD_FORM) else 0)
    else:
      formCount = initialCount

    self.forms = []

    for i in range(0, formCount):
      self.forms.append(FormClass(data=data, prefix='{}-{}'.format(prefix, i)))

    if self.managementForm.cleaned_data.get(ADD_FORM):
      self.forms.append(FormClass(prefix='{}-{}'.format(prefix, formCount)))
    elif self.managementForm.cleaned_data.get(DELETE_FORM):
      self.forms.pop()

  def __iter__(self):
    for f in self.forms:
      yield f

  def is_valid(self):
    return all(f.is_valid() for f in self.forms)

  @property
  def cleaned_data(self):
    return [f.cleaned_data for f in self.forms if f.is_valid()]

class ModelFormSet(FormSet):

  def save(commmit=True):
    return [form.save(commit=commit) for form in self.forms if form.is_valid()]
