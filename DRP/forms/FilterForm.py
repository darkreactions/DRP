from django import forms
from operator import or_, and_

class FilterForm(forms.Form):

  def is_empty(self):
    '''This method fails for foreignkey and many to many fields, for whuich you will need to implement your own methods'''
    if hasattr(self, 'checkFields'):
      keys = self.checkFields
    else:
      keys = self.fields.keys()
    
    return not any(self.cleaned_data.get(key) not in self.fields[key].empty_values for key in keys)

  def fetch(self):
    pass

class FilterFormSet(forms.formsets.BaseFormSet):

  def __init__(self, operator=or_, *args, **kwargs):
    self.extra = 1
    self.can_delete=False
    self.can_order=False
    self.max_num=10000
    self.validate_max=False
    self.absolute_max=10000
    if operator in (or_, and_):
      self._operator = operator
    else:
      raise TypeError('Invalid logical query operator selected, please use operator.or_ or operator.and_')
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
    if hasattr(self, 'model'):
      if self._operator == or_:
        qs = self.model.objects.none()
      elif self._operator == and_:
        qs = self.model.objects.all()
      else:
        raise TypeError('Invalid logical query operator provided') 
    else:
      qs = self.forms[0].fetch()
    for form in self:
      if not form.is_empty():
        qs=self._operator(qs, form.fetch())
    return qs

def filterFormSetFactory(formClass, modelClass):

  class FilterFormSetDerivative(FilterFormSet):

    form = formClass
    model=modelClass

  return FilterFormSetDerivative
