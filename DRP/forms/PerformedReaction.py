'''A module containing code pertaining to forms for the reaction classes'''
from django import forms
from DRP.models import PerformedReaction, RecommendedReaction
from DRP.models import DataSetRelation
from django.contrib.auth.models import User
from django.forms.widgets import HiddenInput

class PerformedRxnAdminForm(forms.ModelForm):
  '''A form for the performed reactions in the Django admin, the only things we want editing from
  there are whether a Reaction is invalid or if it should have the legacy recommended flag set'''

  class Meta:
    fields=('valid', 'legacyRecommendedFlag')
    model=PerformedReaction

class PerformedRxnForm(forms.ModelForm):
  '''A form for creating performed reaction instances in teh databases'''
 
  class Meta:
    fields=('reference', 'performedBy', 'performedDateTime', 'notes', 'labGroup', 'recommendation', 'public', 'duplicateOf', 'valid') 
    #TODO: add tooltip to validity field for initial reaction creation step.
    model=PerformedReaction

  def __init__(self, user, *args, **kwargs):
    '''Overridden __init__ method; requires the user as the first argument so that choice of lab group etc can be validated, as well as to
    track who enters what'''
    super(PerformedRxnForm, self).__init__(*args, **kwargs)
    self.user = user
    labGroups = user.labgroup_set.all()
    self.fields['labGroup'].queryset = labGroups
    self.fields['recommendation'].queryset = RecommendedReaction.objects.filter(labGroup__in=labGroups)
    self.fields['recommendation'].widget = forms.HiddenInput()
    self.fields['duplicateOf'].queryset = PerformedReaction.objects.filter(labGroup__in=labGroups)|PerformedReaction.objects.filter(public=True)
    self.fields['performedBy'].queryset = User.objects.filter(labgroup__in=labGroups)
    self.fields['valid'].initial=False;
    self.fields['valid'].help_text = "This should be left unchecked until descriptor data has been added (after save has been pressed)"
    if labGroups.exists():
      self.fields['labGroup'].empty_label = None 
  
  def clean(self):
    self.cleaned_data = super(PerformedRxnForm, self).clean()
    if 'performedBy' in self.cleaned_data:
        if not self.fields['labGroup'] in (self.user.labgroup_set.all()| self.cleaned_data['performedBy'].labgroup_set.all()):
            raise ValidationError('The selected labGroup does not contain both the inputting and experimental user.', 'invalid_lg')
    return self.cleaned_data
    

  def save(self, commit=True, *args, **kwargs):
    '''Overriden save method automates addition of user that created this instance'''
    rxn = super(PerformedRxnForm, self).save(commit=False)
    rxn.user = self.user
    if commit:
      rxn.save()
    return rxn 

class PerformedRxnDeleteForm(forms.ModelForm):

  class Meta:
    fields=('id',)
    model=PerformedReaction

  def __init__(self, user, *args, **kwargs):
    super(PerformedRxnDeleteForm, self).__init__(*args, **kwargs)
    self.fields['id'] = forms.ModelChoiceField(queryset=PerformedReaction.objects.filter(labGroup__in=user.labgroup_set.all()), initial=self.instance.pk, widget=HiddenInput)

  def clean_id(self):
    if DataSetRelation.objects.filter(reaction__pk=self.cleaned_data['id']).exists():
      raise ValidationError("This reaction is protected from deletion because it is used in one or more reactions or recommendations.")
    return self.cleaned_data['id'] 

  def save(self):
    self.cleaned_data['id'].delete()
    return self.cleaned_data['id']

class PerformedRxnInvalidateForm(forms.ModelForm):

  class Meta:
    fields=('id',)
    model=PerformedReaction

  def __init__(self, user, *args, **kwargs):
    super(PerformedRxnDeleteForm, self).__init__(*args, **kwargs)
    self.fields['id'] = forms.ModelChoiceField(queryset=PerformedReaction.objects.filter(labGroup__in=user.labgroup_set.all()), initial=self.instance.pk, widget=HiddenInput)

  def save(self):
    self.cleaned_data['id'].valid = False
    self.cleaned_data['id'].save()
    return self.cleaned_data['id']
