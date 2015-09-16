'''A module containing code pertaining to forms for the reaction classes'''
from django import forms
from DRP.models import PerformedReaction, RecommendedReaction

class PerformedRxnAdminForm(forms.ModelForm):
  '''A form for the performed reactions in the Django admin, the only things we want editing from
  there are whether a Reaction is invalid or if it should have the legacy recommended flag set'''

  class Meta:
    fields=('valid', 'legacyRecommendedFlag')
    model=PerformedReaction

class PerformedRxnForm(forms.ModelForm):
  '''A form for creating performed reaction instances in teh databases'''
 
  class Meta:
    fields=('reference', 'performedDateTime', 'notes', 'labGroup', 'recommendation', 'public', 'duplicateOf' ) 
    model=PerformedReaction

  def __init__(self, user, *args, **kwargs):
    '''Overridden __init__ method; requires the user as teh first argument so that choice of lab group etc can be validated, as well as to
    track who enters what'''
    super(PerformedRxnForm, self).__init__(*args, **kwargs)
    self.user = user
    labGroups = user.labgroup_set.all()
    self.fields['labGroup'].queryset = labGroups
    self.fields['recommendation'].queryset = RecommendedReaction.objects.filter(labGroup__in=labGroups)
    self.fields['duplicateOf'].queryset = PerformedReaction.objects.filter(labGroup__in=labGroups)|PerformedReaction.objects.filter(public=True)
    if labGroups.exists():
      self.fields['labGroup'].empty_label = None 

  def save(self, commit=True, *args, **kwargs):
    '''Overriden save method automates addition of user that created this instance'''
    rxn = super(PerformedRxnForm, self).save(commit=False)
    rxn.user = self.user
    if commit:
      rxn.save()
    return rxn 
