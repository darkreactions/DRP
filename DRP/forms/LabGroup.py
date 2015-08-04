'''A module containing form classes pertaining to the management of the LabGroup class'''
import django.forms as forms
from django.contrib.auth.hashers import make_password
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
      raise ValidationError('Access Code required')
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

