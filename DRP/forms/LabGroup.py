'''A module containing form classes pertaining to the management of the LabGroup class'''
import django.forms as forms
from django.contrib.auth.hashers import make_password, check_password
from django.conf import settings
from DRP.models import LabGroup
from django.core.exceptions import ValidationError

class LabGroupForm(forms.ModelForm):
  '''This class is for use in the Django admin for creating lab groups.
  Has a Meta class setting the relevant model to LabGroup (defined in the Models subpackage of DRP)

  The class implements one new method, clean_accessCode, and overrides the save method.
  '''

  accessCode=forms.CharField(label='Access Code', widget=forms.PasswordInput, required=False, help_text='''The access code cannot be displayed due to security reasons.
  Entering a new access code here will change the access code. If nothing is entered, it will remain the same.''')
  '''This field is included manually in the form rather than using a ModelForm simple conversion because the accessCode is stored
  as a hash'''

  class Meta:
    model = LabGroup
    fields = ("title", "address", "email", 'accessCode', 'users')

  def clean_accessCode(self):
    '''This method permits the use of old-style LabGroups by either saving a new access code or
    converting the legacy access code (previously stored as a plaintext string) before erasing it.
    '''
    if self.instance.legacy_access_code == '' and self.instance.access_code == '' and self.cleaned_data['accessCode'] == '':
      raise ValidationError('Access Code required', code='no_code')
    elif self.instance.access_code == '' and self.cleaned_data['accessCode'] == '':
      return self.instance.legacy_access_code
    else:
      return self.cleaned_data['accessCode']
    

  def save(self, commit=True):
    '''Saves an instance of the LabGroup, hashing the access code for storage.'''
    labGroup = super(LabGroupForm, self).save(commit=False)
    labGroup.access_code = make_password(self.cleaned_data['accessCode'], settings.LAB_GROUP_HASH_SALT)
    labGroup.legacy_access_code = ''
    if commit:
      labGroup.save()
    return labGroup

class LabGroupJoiningForm(forms.Form):
  '''This class is to validate a user to join a lab group using their supplied access code'''

  labGroup = forms.ModelChoiceField(label='Lab Group', queryset=LabGroup.objects.all())
  accessCode = forms.CharField(label='Access Code', widget=forms.PasswordInput)

  def clean(self):
    super(LabGroupJoiningForm, self).clean()
    if check_password(self.cleaned_data['accessCode'], self.cleaned_data['labGroup'].access_code):
      return self.cleaned_data
    elif self.cleaned_data['accessCode'] == self.cleaned_data['labGroup'].legacy_access_code:
      return self.cleaned_data
    else:
      raise ValidationError('Invalid Access Code', code='invalid_access')

class LabGroupLeavingForm(forms.Form):
  '''This class is a form to allow a user to leave a research group'''

  def __init__(self, user, labGroup = None, *args, **kwargs):
    super(LabGroupLeavingForm, self).__init__(*args, **kwargs)
    self.fields['labGroup'] = forms.ModelChoiceField(queryset=user.labgroup_set.all(), widget=forms.HiddenInput, initial=labGroup)

class LabGroupSelectionForm(forms.Form):
  '''This class is to validate a user to select a group in order to view the compound lists and reaction lists'''
    
  def __init__(self, user, *args, **kwargs):
    super(LabGroupSelectionForm, self).__init__(*args, **kwargs)
    self.fields['labGroup'] = forms.ModelChoiceField(label='Viewing as Lab Group', queryset=user.labgroup_set.all())
