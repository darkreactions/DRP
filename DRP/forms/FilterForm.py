
from django import forms
import abc

class FilterForm(forms.Form):

  __metaclass__ = abc.ABCMeta

  def is_empty(self):
    if hasattr(self, 'checkFields'):
      keys = self.checkFields
    else:
      keys = self.fields.keys()

    return not any(self.cleaned_data.get(key) not in self.fields[key].empty_values for key in keys)

  @abc.abstractmethod
  def fetch(self):
    pass

def FilterFormSet(forms.formsets.BaseFormSet):

  def __init__(self, *args, **kwargs):
    self.extra = 1
    self.can_delete=False
    self.can_order=False
    self.max_num=False
    self.validate_max=False
    self.absolute_max=10000
    super(FilterFormSet, self).__init__(*args, **kwargs)
    
  @property
  def cleaned_data(self):
    initData = []
    for form in self:
      if form.is_valid() and not form.is_empty():
        initData.append(form.cleaned_data)
    return initData 

  def is_empty(self):
    return all(form.is_empty() for form in self)

  def fetch(self):
    qs = self.forms[0].fetch()
    for form in self:
      if not form.is_empty():
        qs |= form.fetch()
    return qs

def filterFormSetFactory(formClass):

  class FilterFormSetDerivative(FilterFormSet):

    form = formClass

  return FilterFormSetDerivative
