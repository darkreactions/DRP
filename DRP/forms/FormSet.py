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
  '''This class acts as a management form for the FormSet class, providing information for the javascript side of things as well as actually functioning in
  line with HTTP standards (which is the half where django's inbuilt formsets seem to not work at present).'''

  def __init__(self, maxForms, formSetPrefix, initialCount, canAdd, canDelete, prefix='', data=None, *args, **kwargs):
    '''Initialiser takes the following arguments:

      maxForms is the maximum number of forms which can be found in this formset
      formSetPrefix is the prefix for the formSet which this class manages, not to be confused with prefix'
      initialCount is the initial number of forms in the managed formSet
      canAdd indicates that additional forms can be added to the managed formset and adds an addition button to this form
      canDelete indicates that forms can be removed from the managed formset and adds a deletion button to this form when there are more than
        the minimum number of forms present
      prefix is the prefix that is assigned to the management form by the managed formset
      data is the form data which has been passed to the formset (normally request.POST)
      *args, **kwargs; other nonpertinent arguments which are handed on to the forms.Form constructor for completeness
    '''
    super(FormSetManagerForm, self).__init__(prefix=prefix+'-'+formSetPrefix, data=copy.copy(data), *args, **kwargs)
    self.fields[TOTAL_FORMS] = forms.IntegerField(min_value=0, widget=forms.widgets.HiddenInput, initial=initialCount if canAdd else initialCount+1) #useful for javascript
    self.fields[MAX_FORMS] = forms.IntegerField(min_value=0, max_value=maxForms, initial=maxForms, widget=forms.widgets.HiddenInput) #useful for javascript
    self.fields[PREFIX] = forms.CharField(initial=formSetPrefix, widget=forms.widgets.HiddenInput, required=False) #useful for javascript
    if canAdd:
      self.canAdd = canAdd
      if data is None:
        self.fields[ADD_FORM] = forms.BooleanField(label=None, widget=SubmitButtonWidget('Add one'), required=False)
      else:
        self.fields[ADD_FORM] = forms.BooleanField(label=None, widget=SubmitButtonWidget('Add another'), required=False)
        
    if canDelete:
      self.canDelete = canDelete
      self.fields[DELETE] = forms.BooleanField(label=None, widget=SubmitButtonWidget('Remove one'), required=False)

  def clean(self):
    '''Overridden clean method checks that the number of forms is correct, and alteres the internal data
    to ensure that it matches the number of forms after any additions have been made'''
    cleaned_data = super(FormSetManagerForm, self).clean()
    if cleaned_data.get(TOTAL_FORMS) > cleaned_data.get(MAX_FORMS):
      raise ValidationError('No more may be added.', code='too_many_forms')
    if cleaned_data.get(ADD_FORM):
      cleaned_data[TOTAL_FORMS] +=1
    if cleaned_data.get(DELETE):
      cleaned_data[TOTAL_FORMS] -=1
    self.data[self.add_prefix(TOTAL_FORMS)] = str(self.cleaned_data.get(TOTAL_FORMS)) #YUK.
    if cleaned_data.get(TOTAL_FORMS) < 1 and (DELETE in self.fields): #this minimum form value should be made a variable at some poitn
      del self.fields[DELETE]
    if cleaned_data.get(TOTAL_FORMS) < 1 and (ADD_FORM in self.fields):
      self.fields[ADD_FORM].widget.value = 'Add one' 
    self.cleaned_data = cleaned_data
    return cleaned_data

  @property
  def submittedForms(self):
    '''Returns the correct number of forms that were actually submitted with data on this occasion'''
    if self.cleaned_data.get(ADD_FORM):
      return self.cleaned_data.get(TOTAL_FORMS)-1
    else:
      return self.cleaned_data.get(TOTAL_FORMS)
      
  def as_ul(self):
    '''Overrides the html output for this form since it only has submit buttons which don't have labels.
    This is horribly hacky but django doesn't actually offer a better way to do this at present'''
    return self._html_output(
      normal_row = '<li%(html_class_attr)s>%(errors)s%(field)s%(help_text)s</li>',
      error_row = '<li>%s</li>',
      row_ender = "</li>",
      help_text_html = '<span class="helptext">%s</span>',
      errors_on_separate_row=False
    )

class FormSet(object):
  '''A drop in replacement for django formsets which more accurately supports HTTP standard way of doing things where javascript cannot be
  depended upon'''

  def __init__(self, formClass, data=None, prefix='', initial=None, maxForms=HARD_MAX_FORMS, canDelete=False, canAdd=True):
    '''initialiser accepts the following arguments:

      formClass is the class of form that is being accepted
    '''
    self.formClass = formClass
    self.prefix=prefix
    self.initial = initial
    self.canAdd = canAdd
    initialCount = self._initialCount()
    self.managementForm = FormSetManagerForm(maxForms, prefix, initialCount, canAdd, canDelete, prefix='{}-manager'.format(prefix), data=data)
    if self.managementForm.is_bound and self.managementForm.is_valid():
        formCount = self.managementForm.submittedForms
    else:
      formCount = initialCount

    self.data = data
    self._createForms(formCount)

    if self.managementForm.is_valid():
      if self.managementForm.cleaned_data.get(ADD_FORM):
        self.forms.append(formClass(prefix='{}-{}'.format(prefix, formCount)))

  def _createForms(self, formCount):
    self.forms = []
    for i in range(0, formCount):
      self.forms.append(self.formClass(data=self.data, prefix='{}-{}'.format(self.prefix, i)))

  def _initialCount(self):
    if self.initial is None:
      if self.canAdd:
        initialCount = 0
      else:
        initalCount = 1
    else:
      initialCount = len(self.initial)
    return initialCount

  def __iter__(self):
    for f in self.forms:
      yield f

  def __len__(self):
    return len(self.forms)

  def is_valid(self):
    '''Checks management form and all forms are valid and returns false otherwise'''
    return self.managementForm.is_valid() and all(f.is_valid() for f in self.forms)

  @property
  def cleaned_data(self):
    '''returns cleaned data for all valid forms'''
    return [f.cleaned_data for f in self.forms if f.is_valid()]

class ModelFormSet(FormSet):
  '''A modified formset to deal with django models'''

  def __init__(self, modelClass, formClass=None, fields=None, instances=None, *args, **kwargs):
    ''' Overridden initialiser accepts the following arguments: 

      modelClass is the django model class that this formset is for
      formClass is the class of form to use if there is one already available,
      if not, then one is automagically constructed.
      fields is an iterable containing the names of the fields which should
      be present. if a formClass has been supplied the fields there are
      overridden.
      *args and **kwargs are the arguments that are suppled to the FormSet class'''
    self.formClass = formClass
    self.instances = instances
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

  def _createForms(self, formCount):
    self.forms = []
    for i in range(0, formCount):
      instance = self.instances[i] if i in range(0, self.instances.count()) else None 
      self.forms.append(self.formClass(data=self.data, instance=instance, prefix='{}-{}'.format(self.prefix, i)))

  def _initialCount(self):
    if self.instances is None:
      if self.canAdd:
        initialCount = 0
      else:
        initalCount = 1
    else:
      initialCount = len(self.instances)
    return initialCount

  def save(self, commit=True):
    '''Returns the objects created by all valid forms in this formset'''
    return [form.save(commit=commit) for form in self.forms if form.is_valid()]
