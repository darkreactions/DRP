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
 
  class Meta:
    fields=('notes', 'labGroup', 'reference', 'recommendation', 'public', 'duplicateOf') 
    model=PerformedReaction

  def __init__(self, user, *args, **kwargs):
    super(PerformedRxnForm, self).__init(*args, **kwargs)
    self.user = user
    labGroups = user.labgroup_set.all()
    self.fields['labGroup'].queryset = labGroups
    self.fields['recommendation'].queryset = RecommendedReaction.objects.filter(labGroup__in=labGroups)
    self.fields['duplicateOf'].queryset = PerformedReactions.objects.filter(labGroup__in=labGroups)|PerformedReactions.objects.filter(public=True)
    if labGroups.exists():
      self.fields['labGroup'].empty_label = None 

  def save(self, commit=True, *args, **kwargs):
    rxn = super(PerformedRxnForm, self).save(commit=False)
    rxn.user = self.user
    if commit:
      rxn.save()
    return rxn 
